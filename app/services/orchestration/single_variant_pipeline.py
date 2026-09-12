from app.models.annotated_variant import AnnotatedVariant
from app.services.annotation.variant_annotator import annotate_variant
from app.services.normalization.variant_normalizer import normalize_variant


def run_single_variant_pipeline(submitted_variant: str) -> AnnotatedVariant:
    normalized_variant = normalize_variant(submitted_variant)
    return annotate_variant(normalized_variant)
