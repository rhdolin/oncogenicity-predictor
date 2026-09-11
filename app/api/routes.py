from fastapi import APIRouter, HTTPException, Query

from app.models.normalized_variant import NormalizedVariant, NormalizedVariantBatchResponse
from app.models.requests import BatchRequest
from app.services.normalization.variant_normalizer import VariantNormalizationError
from app.services.orchestration.single_variant_pipeline import run_single_variant_pipeline


router = APIRouter()


def _normalize_variant_or_raise(submitted_variant: str) -> NormalizedVariant:
    try:
        return run_single_variant_pipeline(submitted_variant)
    except VariantNormalizationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get(
    "/predictOncogenicity",
    response_model=NormalizedVariant,
    summary="Normalize a single variant",
    description=(
        "Accepts one submitted variant in HGVS format, sends it to the ClinGen "
        "Allele Registry for normalization, and returns the internal "
        "NormalizedVariant JSON shape."
    ),
)
def predict_single(
    variant: str = Query(
        description="Variant to normalize. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
) -> NormalizedVariant:
    return _normalize_variant_or_raise(variant)


@router.post(
    "/predictOncogenicity",
    response_model=NormalizedVariantBatchResponse,
    summary="Normalize a batch of variants",
    description=(
        "Accepts a list of submitted variants in HGVS format, normalizes each one "
        "through ClinGen, and returns a list of internal NormalizedVariant objects."
    ),
)
def predict_batch(request: BatchRequest) -> NormalizedVariantBatchResponse:
    normalized_variants = [
        _normalize_variant_or_raise(variant)
        for variant in request.variants
    ]
    return NormalizedVariantBatchResponse(normalized_variants=normalized_variants)
