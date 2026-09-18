"""Run local reference-set evaluation and write prediction and metric CSVs.

This script keeps the evaluation workflow local to the repository. It reads the
selected reference set, runs the existing normalization/annotation/scoring
pipeline for each variant, writes a flattened FHIR-based comparison table, and
generates separate category-based and criteria-based metric summaries.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys


EVALUATION_DIR = Path(__file__).resolve().parent
REPO_ROOT = EVALUATION_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.models.prediction import OncogenicityObservation, OncogenicityPredictionSummary  # noqa: E402
from app.services.fhir import build_oncogenicity_observation  # noqa: E402
from app.services.normalization.variant_normalizer import VariantNormalizationError  # noqa: E402
from app.services.orchestration.single_variant_pipeline import (  # noqa: E402
    run_single_variant_evidence_summary_pipeline,
)


REFERENCE_SET_PATH = EVALUATION_DIR / "variantLists" / "variantList.csv"
OUTPUT_DIR = EVALUATION_DIR / "output"
PREDICTIONS_OUTPUT_PATH = OUTPUT_DIR / "oncogenicityPredictions.csv"
CATEGORY_METRICS_OUTPUT_PATH = OUTPUT_DIR / "metrics-category-based.csv"
SCORE_METRICS_OUTPUT_PATH = OUTPUT_DIR / "metrics-score-based.csv"

REFERENCE_FIELDNAMES = [
    "variant",
    "gene",
    "source",
    "classification",
    "points",
    "criteria",
    "referenceDetailLevel",
    "comments",
]
PREDICTION_FIELDNAMES = [
    "variant",
    "gene",
    "source",
    "referenceDetailLevel",
    "referenceClassification",
    "referenceScore",
    "referenceCriteria",
    "predictedClassification",
    "predictedScore",
    "predictedCriteria",
    "predictionUnavailable",
    "evidenceSummary",
]
SWIMLANE_ORDER = ["population", "computational", "hotspots", "predictive", "om1", "op2", "functional"]
CRITERIA_MAP = {
    "SBVS1": ("population", -8),
    "SBS1": ("population", -4),
    "OP4": ("population", 1),
    "SBP1": ("computational", -1),
    "OP1": ("computational", 1),
    "OP3": ("hotspots", 1),
    "OM3": ("hotspots", 2),
    "OS3": ("hotspots", 4),
    "SBP2": ("predictive", -1),
    "OM2": ("predictive", 2),
    "OM4": ("predictive", 2),
    "OS1": ("predictive", 4),
    "OVS1": ("predictive", 8),
    "OM1": ("om1", 2),
    "OP2": ("op2", 1),
    "SBS2": ("functional", -4),
    "OS2-2": ("functional", 2),
    "OS2": ("functional", 4),
}
CLASSIFICATION_NORMALIZATION_MAP = {
    "benign": "Benign",
    "likely benign": "Likely Benign",
    "vus": "VUS",
    "likely oncogenic": "Likely Oncogenic",
    "oncogenic": "Oncogenic",
}
REFERENCE_DETAIL_LEVELS = {"full", "category_only"}


def load_reference_rows(reference_set_path: Path = REFERENCE_SET_PATH) -> list[dict[str, str]]:
    with reference_set_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing_columns = [column for column in REFERENCE_FIELDNAMES if column not in (reader.fieldnames or [])]
        if missing_columns:
            raise ValueError(f"Reference set is missing required columns: {missing_columns}")

        rows: list[dict[str, str]] = []
        for raw_row in reader:
            variant = (raw_row.get("variant") or "").strip()
            if not variant:
                continue

            rows.append(
                {
                    field_name: (raw_row.get(field_name) or "").strip()
                    for field_name in REFERENCE_FIELDNAMES
                }
            )
    return rows


def normalize_reference_detail_level(raw_value: str) -> str:
    normalized_value = raw_value.strip().lower()
    if normalized_value not in REFERENCE_DETAIL_LEVELS:
        raise ValueError(
            "referenceDetailLevel must be one of: category_only, full"
        )
    return normalized_value


def parse_reference_score(raw_value: str, reference_detail_level: str) -> int | None:
    stripped_value = raw_value.strip()
    if not stripped_value:
        if reference_detail_level == "category_only":
            return None
        raise ValueError("Reference score is required for full-detail rows.")
    if reference_detail_level == "category_only":
        raise ValueError("Category-only rows must leave points blank.")
    return int(stripped_value)


def parse_criteria_list(raw_value: str, reference_detail_level: str) -> list[str]:
    stripped_value = raw_value.strip()
    if not stripped_value:
        if reference_detail_level == "category_only":
            return []
        raise ValueError("Reference criteria are required for full-detail rows.")
    if reference_detail_level == "category_only":
        raise ValueError("Category-only rows must leave criteria blank.")

    parsed_value = json.loads(stripped_value)

    if not isinstance(parsed_value, list):
        raise ValueError(f"Expected criteria list, got: {type(parsed_value).__name__}")

    criteria: list[str] = []
    for item in parsed_value:
        if not isinstance(item, str):
            continue
        normalized_item = item.strip()
        if normalized_item:
            criteria.append(normalized_item)
    return criteria


def build_reference_fields(reference_row: dict[str, str]) -> dict[str, object]:
    reference_detail_level = normalize_reference_detail_level(
        reference_row["referenceDetailLevel"]
    )
    return {
        "referenceDetailLevel": reference_detail_level,
        "referenceClassification": normalize_classification(reference_row["classification"]),
        "referenceScore": parse_reference_score(
            reference_row["points"],
            reference_detail_level,
        ),
        "referenceCriteria": parse_criteria_list(
            reference_row["criteria"],
            reference_detail_level,
        ),
    }


def normalize_classification(raw_value: str) -> str:
    normalized_key = raw_value.strip().lower()
    if not normalized_key:
        return ""
    return CLASSIFICATION_NORMALIZATION_MAP.get(normalized_key, raw_value.strip())


def build_prediction_row(reference_row: dict[str, str]) -> dict[str, object]:
    variant = reference_row["variant"]
    reference_fields = build_reference_fields(reference_row)
    try:
        summary = run_single_variant_evidence_summary_pipeline(variant)
    except VariantNormalizationError as exc:
        return {
            "variant": variant,
            "gene": reference_row["gene"],
            "source": reference_row["source"],
            **reference_fields,
            "predictedClassification": "",
            "predictedScore": "",
            "predictedCriteria": [],
            "predictionUnavailable": True,
            "evidenceSummary": {"detail": str(exc)},
        }

    observation = build_oncogenicity_observation(summary, variant)
    return flatten_prediction_row(reference_row, summary, observation)


def flatten_prediction_row(
    reference_row: dict[str, str],
    summary: OncogenicityPredictionSummary,
    observation: OncogenicityObservation,
) -> dict[str, object]:
    reference_fields = build_reference_fields(reference_row)
    predicted_classification = ""
    if observation.interpretation:
        predicted_classification = (
            observation.interpretation[0].coding[0].code
            or observation.interpretation[0].text
            or ""
        )

    predicted_score: int | str = ""
    if observation.valueInteger is not None:
        predicted_score = observation.valueInteger

    predicted_criteria = extract_predicted_criteria(observation)

    return {
        "variant": reference_row["variant"],
        "gene": reference_row["gene"],
        "source": reference_row["source"],
        **reference_fields,
        "predictedClassification": predicted_classification,
        "predictedScore": predicted_score,
        "predictedCriteria": predicted_criteria,
        "predictionUnavailable": observation.dataAbsentReason is not None,
        "evidenceSummary": summary.model_dump(mode="json"),
    }


def extract_predicted_criteria(observation: OncogenicityObservation) -> list[str]:
    criteria: list[str] = []
    for component in observation.component:
        if not component.interpretation:
            continue
        concept = component.interpretation[0]
        if concept.coding and concept.coding[0].code:
            criteria.append(concept.coding[0].code)
            continue
        if concept.text:
            criteria.append(concept.text)
    return criteria


def serialize_prediction_row(row: dict[str, object]) -> dict[str, object]:
    return {
        "variant": row["variant"],
        "gene": row["gene"],
        "source": row["source"],
        "referenceDetailLevel": row["referenceDetailLevel"],
        "referenceClassification": row["referenceClassification"],
        "referenceScore": "" if row["referenceScore"] is None else row["referenceScore"],
        "referenceCriteria": json.dumps(row["referenceCriteria"]),
        "predictedClassification": row["predictedClassification"],
        "predictedScore": row["predictedScore"],
        "predictedCriteria": json.dumps(row["predictedCriteria"]),
        "predictionUnavailable": json.dumps(row["predictionUnavailable"]),
        "evidenceSummary": json.dumps(row["evidenceSummary"]),
    }


def write_predictions_csv(rows: list[dict[str, object]], output_path: Path = PREDICTIONS_OUTPUT_PATH) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREDICTION_FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(serialize_prediction_row(row))


def generate_category_metrics(rows: list[dict[str, object]]) -> tuple[dict[str, float | int], list[list[int]]]:
    available_rows = [row for row in rows if not row["predictionUnavailable"]]
    num_input_rows = len(rows)
    num_prediction_unavailable = num_input_rows - len(available_rows)
    num_category_tested = len(available_rows)
    score_testable_rows = [
        row for row in available_rows if row["referenceDetailLevel"] == "full"
    ]
    num_score_tested = len(score_testable_rows)

    benign_to_benign = benign_to_vus = benign_to_oncogenic = 0
    vus_to_benign = vus_to_vus = vus_to_oncogenic = 0
    oncogenic_to_benign = oncogenic_to_vus = oncogenic_to_oncogenic = 0
    exact_score_agreement = score_within_2 = score_greater_than_2 = 0

    for row in available_rows:
        reference_bucket = classification_bucket(str(row["referenceClassification"]))
        predicted_bucket = classification_bucket(str(row["predictedClassification"]))

        if reference_bucket == "Benign":
            if predicted_bucket == "Benign":
                benign_to_benign += 1
            elif predicted_bucket == "VUS":
                benign_to_vus += 1
            else:
                benign_to_oncogenic += 1
        elif reference_bucket == "VUS":
            if predicted_bucket == "Benign":
                vus_to_benign += 1
            elif predicted_bucket == "VUS":
                vus_to_vus += 1
            else:
                vus_to_oncogenic += 1
        else:
            if predicted_bucket == "Benign":
                oncogenic_to_benign += 1
            elif predicted_bucket == "VUS":
                oncogenic_to_vus += 1
            else:
                oncogenic_to_oncogenic += 1

    for row in score_testable_rows:
        predicted_score = int(row["predictedScore"])
        reference_score = int(row["referenceScore"])
        absolute_delta = abs(reference_score - predicted_score)
        if absolute_delta == 0:
            exact_score_agreement += 1
        if absolute_delta <= 2:
            score_within_2 += 1
        else:
            score_greater_than_2 += 1

    adjacent_disagreements = benign_to_vus + vus_to_benign + vus_to_oncogenic + oncogenic_to_vus
    extreme_disagreements = benign_to_oncogenic + oncogenic_to_benign
    diagonal_agreement = benign_to_benign + vus_to_vus + oncogenic_to_oncogenic

    row_benign = benign_to_benign + benign_to_vus + benign_to_oncogenic
    row_vus = vus_to_benign + vus_to_vus + vus_to_oncogenic
    row_oncogenic = oncogenic_to_benign + oncogenic_to_vus + oncogenic_to_oncogenic
    column_benign = benign_to_benign + vus_to_benign + oncogenic_to_benign
    column_vus = benign_to_vus + vus_to_vus + oncogenic_to_vus
    column_oncogenic = benign_to_oncogenic + vus_to_oncogenic + oncogenic_to_oncogenic

    observed_agreement = safe_rate(diagonal_agreement, num_category_tested)
    weighted_agreement_linear = safe_rate(
        diagonal_agreement + 0.5 * adjacent_disagreements,
        num_category_tested,
    )
    weighted_agreement_quadratic = safe_rate(
        diagonal_agreement + 0.75 * adjacent_disagreements,
        num_category_tested,
    )

    expected_weighted_linear = safe_rate(
        (
            1.0 * row_benign * column_benign
            + 0.5 * row_benign * column_vus
            + 0.5 * row_vus * column_benign
            + 1.0 * row_vus * column_vus
            + 0.5 * row_vus * column_oncogenic
            + 0.5 * row_oncogenic * column_vus
            + 1.0 * row_oncogenic * column_oncogenic
        ),
        num_category_tested * num_category_tested,
    )
    expected_weighted_quadratic = safe_rate(
        (
            1.0 * row_benign * column_benign
            + 0.75 * row_benign * column_vus
            + 0.75 * row_vus * column_benign
            + 1.0 * row_vus * column_vus
            + 0.75 * row_vus * column_oncogenic
            + 0.75 * row_oncogenic * column_vus
            + 1.0 * row_oncogenic * column_oncogenic
        ),
        num_category_tested * num_category_tested,
    )

    metrics = {
        "numInputRows": num_input_rows,
        "numPredictionUnavailable": num_prediction_unavailable,
        "numCategoryTested": num_category_tested,
        "numScoreTested": num_score_tested,
        "observedAgreement": observed_agreement,
        "weightedAgreementLinear": weighted_agreement_linear,
        "expectedWeightedLinear": expected_weighted_linear,
        "weightedKappaLinear": safe_kappa(weighted_agreement_linear, expected_weighted_linear),
        "weightedAgreementQuadratic": weighted_agreement_quadratic,
        "expectedWeightedQuadratic": expected_weighted_quadratic,
        "weightedKappaQuadratic": safe_kappa(
            weighted_agreement_quadratic,
            expected_weighted_quadratic,
        ),
        "adjacentDisagreements": adjacent_disagreements,
        "extremeDisagreements": extreme_disagreements,
        "diagonalAgreement": diagonal_agreement,
        "exactScoreAgreement": exact_score_agreement,
        "exactScoreAgreementRate": safe_rate(exact_score_agreement, num_score_tested),
        "scoreWithin2": score_within_2,
        "scoreWithin2Rate": safe_rate(score_within_2, num_score_tested),
        "scoreGreaterThan2": score_greater_than_2,
        "scoreGreaterThan2Rate": safe_rate(score_greater_than_2, num_score_tested),
    }

    confusion_matrix = [
        [benign_to_benign, benign_to_vus, benign_to_oncogenic],
        [vus_to_benign, vus_to_vus, vus_to_oncogenic],
        [oncogenic_to_benign, oncogenic_to_vus, oncogenic_to_oncogenic],
    ]
    return metrics, confusion_matrix


def write_category_metrics(
    metrics: dict[str, float | int],
    confusion_matrix: list[list[int]],
    output_path: Path = CATEGORY_METRICS_OUTPUT_PATH,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics.keys()))
        writer.writeheader()
        writer.writerow(metrics)
        handle.write("\n")
        handle.write("Confusion Matrix (Reference -> Predicted)\n")
        writer = csv.writer(handle)
        writer.writerow(["", "Predicted_Benign", "Predicted_VUS", "Predicted_Oncogenic"])
        for label, row in zip(
            ["Reference_Benign", "Reference_VUS", "Reference_Oncogenic"],
            confusion_matrix,
            strict=True,
        ):
            writer.writerow([label, *row])


def generate_score_metrics(
    rows: list[dict[str, object]],
) -> tuple[dict[str, float | int], list[dict[str, float | int | str]]]:
    available_rows = [row for row in rows if not row["predictionUnavailable"]]
    criteria_testable_rows = [
        row for row in available_rows if row["referenceDetailLevel"] == "full"
    ]
    num_input_rows = len(rows)
    num_prediction_unavailable = num_input_rows - len(available_rows)
    num_criteria_tested = len(criteria_testable_rows)

    exact_criteria_set_matches = 0
    swimlane_summary = {
        swimlane: {"same": 0, "within2": 0, "greaterThan2": 0}
        for swimlane in SWIMLANE_ORDER
    }

    for row in criteria_testable_rows:
        reference_criteria = list(row["referenceCriteria"])
        predicted_criteria = list(row["predictedCriteria"])
        if reference_criteria == predicted_criteria:
            exact_criteria_set_matches += 1

        reference_scores = criteria_to_swimlane_scores(reference_criteria)
        predicted_scores = criteria_to_swimlane_scores(predicted_criteria)
        for swimlane in SWIMLANE_ORDER:
            absolute_delta = abs(reference_scores[swimlane] - predicted_scores[swimlane])
            if absolute_delta == 0:
                swimlane_summary[swimlane]["same"] += 1
            elif absolute_delta <= 2:
                swimlane_summary[swimlane]["within2"] += 1
            else:
                swimlane_summary[swimlane]["greaterThan2"] += 1

    overview = {
        "numInputRows": num_input_rows,
        "numPredictionUnavailable": num_prediction_unavailable,
        "numCriteriaTested": num_criteria_tested,
        "exactCriteriaSetMatches": exact_criteria_set_matches,
        "exactCriteriaSetAgreementRate": safe_rate(exact_criteria_set_matches, num_criteria_tested),
    }

    swimlane_rows: list[dict[str, float | int | str]] = []
    for swimlane in SWIMLANE_ORDER:
        same = swimlane_summary[swimlane]["same"]
        within2 = swimlane_summary[swimlane]["within2"]
        greater_than_2 = swimlane_summary[swimlane]["greaterThan2"]
        exact_or_close = same + within2
        swimlane_rows.append(
            {
                "numTested": num_criteria_tested,
                "numCriteriaTested": num_criteria_tested,
                "swimlane": swimlane,
                "same": same,
                "within2": within2,
                "greaterThan2": greater_than_2,
                "sumCheck": same + within2 + greater_than_2,
                "exactAgreementRate": safe_rate(same, num_criteria_tested),
                "acceptableAgreementRate": safe_rate(exact_or_close, num_criteria_tested),
                "largeDeviationRate": safe_rate(greater_than_2, num_criteria_tested),
                "weightedScoreLinear": safe_rate(same + 0.5 * within2, num_criteria_tested),
                "fuzzinessRatio": safe_rate(within2, exact_or_close),
                "failureToExactRatio": safe_rate(greater_than_2, same),
            }
        )

    return overview, swimlane_rows


def criteria_to_swimlane_scores(criteria: list[str]) -> dict[str, int]:
    scores = {swimlane: 0 for swimlane in SWIMLANE_ORDER}
    for criterion in criteria:
        mapping = CRITERIA_MAP.get(criterion)
        if mapping is None:
            continue
        swimlane, score = mapping
        scores[swimlane] = score
    return scores


def write_score_metrics(
    overview: dict[str, float | int],
    swimlane_rows: list[dict[str, float | int | str]],
    output_path: Path = SCORE_METRICS_OUTPUT_PATH,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        overview_writer = csv.DictWriter(handle, fieldnames=list(overview.keys()))
        overview_writer.writeheader()
        overview_writer.writerow(overview)
        handle.write("\n")
        swimlane_writer = csv.DictWriter(handle, fieldnames=list(swimlane_rows[0].keys()))
        swimlane_writer.writeheader()
        for row in swimlane_rows:
            swimlane_writer.writerow(row)


def classification_bucket(classification: str) -> str:
    normalized = normalize_classification(classification)
    if normalized in {"Benign", "Likely Benign"}:
        return "Benign"
    if normalized == "VUS":
        return "VUS"
    if normalized in {"Likely Oncogenic", "Oncogenic"}:
        return "Oncogenic"
    raise ValueError(f"Unsupported classification for metric bucketing: {classification!r}")


def safe_rate(numerator: float | int, denominator: float | int) -> float:
    if denominator == 0:
        return 0.0
    return float(numerator) / float(denominator)


def safe_kappa(observed: float, expected: float) -> float:
    if expected == 1.0:
        return 0.0
    return (observed - expected) / (1.0 - expected)


def build_prediction_rows_with_progress(
    reference_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    total_rows = len(reference_rows)
    prediction_rows: list[dict[str, object]] = []

    print(f"Starting evaluation for {total_rows} variants...", flush=True)
    for index, reference_row in enumerate(reference_rows, start=1):
        variant = reference_row["variant"]
        print(f"[{index}/{total_rows}] Analyzing {variant}", flush=True)
        prediction_row = build_prediction_row(reference_row)
        prediction_rows.append(prediction_row)

        status_label = "prediction unavailable" if prediction_row["predictionUnavailable"] else "complete"
        print(f"[{index}/{total_rows}] Finished {variant} ({status_label})", flush=True)

    return prediction_rows


def main() -> None:
    reference_rows = load_reference_rows()
    prediction_rows = build_prediction_rows_with_progress(reference_rows)

    write_predictions_csv(prediction_rows)

    category_metrics, category_confusion_matrix = generate_category_metrics(prediction_rows)
    write_category_metrics(category_metrics, category_confusion_matrix)

    score_metrics_overview, score_metrics_rows = generate_score_metrics(prediction_rows)
    write_score_metrics(score_metrics_overview, score_metrics_rows)

    unavailable_count = sum(1 for row in prediction_rows if row["predictionUnavailable"])
    print(
        f"Completed evaluation for {len(prediction_rows)} variants; "
        f"{unavailable_count} predictions were unavailable.",
        flush=True,
    )
    print(f"Wrote {PREDICTIONS_OUTPUT_PATH}")
    print(f"Wrote {CATEGORY_METRICS_OUTPUT_PATH}")
    print(f"Wrote {SCORE_METRICS_OUTPUT_PATH}")


if __name__ == "__main__":
    main()
