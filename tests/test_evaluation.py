import json

from app.models.prediction import (
    CodeableConcept,
    Coding,
    EvidenceResult,
    OncogenicityEvidence,
    OncogenicityObservation,
    OncogenicityPredictionSummary,
    ObservationComponent,
)
from evaluation.runEvaluation import (
    flatten_prediction_row,
    generate_category_metrics,
    generate_score_metrics,
    parse_criteria_list,
    parse_reference_score,
    serialize_prediction_row,
)


def _make_summary(overall_score: int | None, overall_classification: str | None) -> OncogenicityPredictionSummary:
    evidence = OncogenicityEvidence(
        population=EvidenceResult(
            score=1,
            evidenceCode="OP4",
            evidenceStatement="Population statement",
            status="applied",
        ),
        computational=EvidenceResult(
            score=1,
            evidenceCode="OP1",
            evidenceStatement="Computational statement",
            status="applied",
        ),
        hotspots=EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="No hotspot evidence",
            status="applied",
        ),
        predictive=EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="No predictive evidence",
            status="applied",
        ),
        om1=EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="OM1 unavailable",
            status="not_available",
            dataAbsentReason="unsupported",
        ),
        op2=EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="OP2 unavailable",
            status="not_available",
            dataAbsentReason="unsupported",
        ),
        functional=EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="No functional evidence",
            status="applied",
        ),
    )
    return OncogenicityPredictionSummary(
        overallScore=overall_score,
        overallClassification=overall_classification,
        predictionStatement=None,
        dataAbsentReason=None,
        oncogenicityEvidence=evidence,
    )


def _make_observation(
    score: int | None,
    classification: str | None,
    criteria: list[str],
    prediction_unavailable: bool = False,
) -> OncogenicityObservation:
    components = [
        ObservationComponent(
            code=CodeableConcept(
                coding=[Coding(system="test", code=f"lane-{index}", display=f"lane-{index}")],
                text=f"lane-{index}",
            ),
            valueInteger=index,
            interpretation=[
                CodeableConcept(coding=[Coding(system="test", code=criterion, display=criterion)], text=criterion)
            ],
        )
        for index, criterion in enumerate(criteria, start=1)
    ]

    interpretation = []
    if classification is not None:
        interpretation = [
            CodeableConcept(
                coding=[Coding(system="test", code=classification, display=classification)],
                text=classification,
            )
        ]

    return OncogenicityObservation(
        issued="2026-01-01T00:00:00Z",
        code=CodeableConcept(coding=[Coding(system="test", code="oncogenicity-prediction")]),
        valueInteger=score,
        interpretation=interpretation,
        dataAbsentReason=(
            CodeableConcept(coding=[Coding(system="test", code="error", display="error")], text="Unavailable")
            if prediction_unavailable
            else None
        ),
        component=components,
    )


def test_flatten_prediction_row_retains_summary_and_criteria() -> None:
    reference_row = {
        "variant": "NM_004119.3:c.2073T>G",
        "gene": "FLT3",
        "source": "paper",
        "classification": "VUS",
        "points": "2",
        "criteria": '["OP4", "OP1"]',
        "referenceDetailLevel": "full",
    }
    summary = _make_summary(overall_score=2, overall_classification="VUS")
    observation = _make_observation(score=2, classification="VUS", criteria=["OP4", "OP1"])

    row = flatten_prediction_row(reference_row, summary, observation)
    serialized_row = serialize_prediction_row(row)

    assert row["predictedCriteria"] == ["OP4", "OP1"]
    assert row["referenceCriteria"] == ["OP4", "OP1"]
    assert row["predictionUnavailable"] is False
    assert row["evidenceSummary"] == summary.model_dump(mode="json")
    assert json.loads(str(serialized_row["referenceCriteria"])) == ["OP4", "OP1"]
    assert json.loads(str(serialized_row["predictedCriteria"])) == ["OP4", "OP1"]


def test_generate_category_metrics_excludes_prediction_unavailable_rows() -> None:
    rows = [
        {
            "referenceDetailLevel": "full",
            "referenceClassification": "Oncogenic",
            "referenceScore": 10,
            "predictedClassification": "Oncogenic",
            "predictedScore": 10,
            "referenceCriteria": ["OS1", "OP4", "OP1"],
            "predictedCriteria": ["OS1", "OP4", "OP1"],
            "predictionUnavailable": False,
        },
        {
            "referenceDetailLevel": "full",
            "referenceClassification": "Benign",
            "referenceScore": -9,
            "predictedClassification": "",
            "predictedScore": "",
            "referenceCriteria": ["SBVS1", "SBS2", "SBP1"],
            "predictedCriteria": [],
            "predictionUnavailable": True,
        },
    ]

    metrics, confusion_matrix = generate_category_metrics(rows)

    assert metrics["numInputRows"] == 2
    assert metrics["numPredictionUnavailable"] == 1
    assert metrics["numCategoryTested"] == 1
    assert metrics["numScoreTested"] == 1
    assert metrics["observedAgreement"] == 1.0
    assert metrics["exactScoreAgreement"] == 1
    assert confusion_matrix == [[0, 0, 0], [0, 0, 0], [0, 0, 1]]


def test_generate_category_metrics_includes_category_only_rows_but_excludes_them_from_score_metrics() -> None:
    rows = [
        {
            "referenceDetailLevel": "category_only",
            "referenceClassification": "Likely Oncogenic",
            "referenceScore": None,
            "predictedClassification": "Likely Oncogenic",
            "predictedScore": 6,
            "referenceCriteria": [],
            "predictedCriteria": ["OS1", "OP4", "OP1"],
            "predictionUnavailable": False,
        },
        {
            "referenceDetailLevel": "full",
            "referenceClassification": "VUS",
            "referenceScore": 2,
            "predictedClassification": "VUS",
            "predictedScore": 2,
            "referenceCriteria": ["OP4", "OP1"],
            "predictedCriteria": ["OP4", "OP1"],
            "predictionUnavailable": False,
        },
    ]

    metrics, confusion_matrix = generate_category_metrics(rows)

    assert metrics["numCategoryTested"] == 2
    assert metrics["numScoreTested"] == 1
    assert metrics["observedAgreement"] == 1.0
    assert metrics["exactScoreAgreement"] == 1
    assert confusion_matrix == [[0, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_generate_score_metrics_tracks_exact_criteria_set_and_swimlane_scores() -> None:
    rows = [
        {
            "referenceDetailLevel": "full",
            "referenceClassification": "VUS",
            "referenceScore": 2,
            "predictedClassification": "VUS",
            "predictedScore": 2,
            "referenceCriteria": ["OP4", "OP1"],
            "predictedCriteria": ["OP4", "OP1"],
            "predictionUnavailable": False,
        },
        {
            "referenceDetailLevel": "full",
            "referenceClassification": "Likely Oncogenic",
            "referenceScore": 6,
            "predictedClassification": "Likely Oncogenic",
            "predictedScore": 6,
            "referenceCriteria": ["OP4", "OS1", "OP1"],
            "predictedCriteria": ["OP4", "OM4", "OP1"],
            "predictionUnavailable": False,
        },
    ]

    overview, swimlane_rows = generate_score_metrics(rows)
    rows_by_swimlane = {row["swimlane"]: row for row in swimlane_rows}

    assert overview["numCriteriaTested"] == 2
    assert overview["exactCriteriaSetMatches"] == 1
    assert overview["exactCriteriaSetAgreementRate"] == 0.5
    assert rows_by_swimlane["population"]["same"] == 2
    assert rows_by_swimlane["computational"]["same"] == 2
    assert rows_by_swimlane["predictive"]["within2"] == 1
    assert rows_by_swimlane["predictive"]["same"] == 1


def test_generate_score_metrics_excludes_rows_without_reference_criteria() -> None:
    rows = [
        {
            "referenceDetailLevel": "category_only",
            "referenceClassification": "Likely Oncogenic",
            "referenceScore": None,
            "predictedClassification": "Likely Oncogenic",
            "predictedScore": 6,
            "referenceCriteria": [],
            "predictedCriteria": ["OS1", "OP4", "OP1"],
            "predictionUnavailable": False,
        },
        {
            "referenceDetailLevel": "full",
            "referenceClassification": "VUS",
            "referenceScore": 2,
            "predictedClassification": "VUS",
            "predictedScore": 2,
            "referenceCriteria": ["OP4", "OP1"],
            "predictedCriteria": ["OP4", "OP1"],
            "predictionUnavailable": False,
        },
    ]

    overview, swimlane_rows = generate_score_metrics(rows)
    rows_by_swimlane = {row["swimlane"]: row for row in swimlane_rows}

    assert overview["numCriteriaTested"] == 1
    assert overview["exactCriteriaSetMatches"] == 1
    assert rows_by_swimlane["population"]["same"] == 1
    assert rows_by_swimlane["computational"]["same"] == 1


def test_category_only_rows_require_blank_points_and_criteria() -> None:
    assert parse_reference_score("", "category_only") is None
    assert parse_criteria_list("", "category_only") == []


def test_full_rows_require_points_and_criteria() -> None:
    assert parse_reference_score("2", "full") == 2
    assert parse_criteria_list('["OP4"]', "full") == ["OP4"]
