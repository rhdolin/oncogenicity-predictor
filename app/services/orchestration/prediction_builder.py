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
from app.services.scoring import calculate_overall_score


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

    return OncogenicityPredictionSummary(
        overallScore=calculate_overall_score(evidence),
        overallClassification=None,
        oncogenicityEvidence=evidence,
    )
