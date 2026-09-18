"""Build the internal evidence summary from an annotated variant.

This module is the boundary between annotation output and scoring output. It
collects pipeline-specific evidence results, computes the overall score, and
returns the intermediate summary used by both debug and FHIR-facing endpoints.
"""

from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import (
    OncogenicityEvidence,
    OncogenicityPredictionSummary,
)
from app.models.tumor_types import TumorType
from app.services.evidence import (
    build_computational_evidence,
    build_functional_evidence,
    build_hotspots_evidence,
    build_om1_evidence,
    build_op2_evidence,
    build_population_evidence,
    build_predictive_evidence,
)
from app.services.scoring import (
    apply_evidence_interaction_rules,
    calculate_overall_score,
    classify_overall_score,
)


def build_prediction_summary_from_annotated_variant(
    annotated_variant: AnnotatedVariant,
    tumor_type: TumorType | None = None,
) -> OncogenicityPredictionSummary:
    """Assemble all pipeline evidence results and compute the current overall score."""
    evidence = OncogenicityEvidence(
        population=build_population_evidence(annotated_variant),
        computational=build_computational_evidence(annotated_variant),
        hotspots=build_hotspots_evidence(annotated_variant),
        predictive=build_predictive_evidence(
            annotated_variant,
            tumor_type=tumor_type,
        ),
        om1=build_om1_evidence(annotated_variant),
        op2=build_op2_evidence(annotated_variant, tumor_type=tumor_type),
        functional=build_functional_evidence(annotated_variant, tumor_type=tumor_type),
    )

    adjusted_evidence = apply_evidence_interaction_rules(evidence)
    if _all_evidence_not_available(adjusted_evidence):
        return _build_unavailable_prediction_summary(annotated_variant, adjusted_evidence)

    overall_score = calculate_overall_score(adjusted_evidence)

    return OncogenicityPredictionSummary(
        overallScore=overall_score,
        overallClassification=classify_overall_score(overall_score),
        predictionStatement="Overall oncogenicity prediction computed from the available evidence.",
        dataAbsentReason=None,
        oncogenicityEvidence=adjusted_evidence,
    )


def _all_evidence_not_available(evidence: OncogenicityEvidence) -> bool:
    return all(
        result.status == "not_available"
        for result in (
            evidence.population,
            evidence.computational,
            evidence.hotspots,
            evidence.predictive,
            evidence.om1,
            evidence.op2,
            evidence.functional,
        )
    )


def _build_unavailable_prediction_summary(
    annotated_variant: AnnotatedVariant,
    evidence: OncogenicityEvidence,
) -> OncogenicityPredictionSummary:
    if annotated_variant.annotationStatus != "complete":
        statement = (
            "Overall oncogenicity prediction could not be determined because variant annotation failed."
        )
        absent_reason = "error"
    else:
        statement = (
            "Overall oncogenicity prediction could not be determined because no evidence pipelines were evaluable."
        )
        absent_reason = "unsupported"

    return OncogenicityPredictionSummary(
        overallScore=None,
        overallClassification=None,
        predictionStatement=statement,
        dataAbsentReason=absent_reason,
        oncogenicityEvidence=evidence,
    )
