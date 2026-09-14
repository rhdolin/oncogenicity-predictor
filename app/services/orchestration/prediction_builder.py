from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import (
    OncogenicityEvidence,
    OncogenicityPredictionSummary,
)
from app.services.evidence import (
    build_computational_evidence,
    build_not_available_evidence,
    build_population_evidence,
)
from app.services.scoring import calculate_overall_score


def build_prediction_summary_from_annotated_variant(
    annotated_variant: AnnotatedVariant,
) -> OncogenicityPredictionSummary:
    evidence = OncogenicityEvidence(
        population=build_population_evidence(annotated_variant),
        computational=build_computational_evidence(annotated_variant),
        hotspots=build_not_available_evidence(
            "Hotspots evidence is not yet implemented in the prediction pipeline."
        ),
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
