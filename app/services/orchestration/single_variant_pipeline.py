from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import OncogenicityObservation
from app.services.annotation.variant_annotator import annotate_variant
from app.services.fhir import build_oncogenicity_observation
from app.services.normalization.variant_normalizer import normalize_variant

from .prediction_builder import build_prediction_summary_from_annotated_variant


def run_single_variant_annotation_pipeline(submitted_variant: str) -> AnnotatedVariant:
    normalized_variant = normalize_variant(submitted_variant)
    return annotate_variant(normalized_variant)


def run_single_variant_pipeline(submitted_variant: str) -> OncogenicityObservation:
    annotated_variant = run_single_variant_annotation_pipeline(submitted_variant)
    summary = build_prediction_summary_from_annotated_variant(annotated_variant)
    return build_oncogenicity_observation(summary, submitted_variant)
