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
from app.services.evidence import (
    build_computational_evidence,
    build_hotspots_evidence,
    build_not_available_evidence,
    build_population_evidence,
)
from app.services.scoring import calculate_overall_score


def build_prediction_summary_from_annotated_variant(
    annotated_variant: AnnotatedVariant,
) -> OncogenicityPredictionSummary:
    """Assemble all pipeline evidence results and compute the current overall score."""
    evidence = OncogenicityEvidence(
        population=build_population_evidence(annotated_variant),
        computational=build_computational_evidence(annotated_variant),
        hotspots=build_hotspots_evidence(annotated_variant),
        predictive=build_not_available_evidence(
            "Predictive evidence is not yet implemented in the prediction pipeline."
        ),
        functional=build_not_available_evidence(
            "Functional evidence is not yet implemented in the prediction pipeline."
        ),
    )

    return OncogenicityPredictionSummary(
        overallScore=calculate_overall_score(evidence),
        overallClassification=None,
        oncogenicityEvidence=evidence,
    )
