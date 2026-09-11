from app.models.normalized_variant import NormalizedVariant
from app.services.normalization.variant_normalizer import normalize_variant


def run_single_variant_pipeline(submitted_variant: str) -> NormalizedVariant:
    return normalize_variant(submitted_variant)
