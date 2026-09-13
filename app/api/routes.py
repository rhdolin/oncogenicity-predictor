from fastapi import APIRouter, HTTPException, Query

from app.models.annotated_variant import AnnotatedVariant, AnnotatedVariantBatchResponse
from app.models.requests import BatchRequest
from app.services.normalization.variant_normalizer import VariantNormalizationError
from app.services.orchestration.single_variant_pipeline import run_single_variant_pipeline


router = APIRouter()


def _annotate_variant_or_raise(submitted_variant: str) -> AnnotatedVariant:
    try:
        return run_single_variant_pipeline(submitted_variant)
    except VariantNormalizationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get(
    "/predictOncogenicity",
    response_model=AnnotatedVariant,
    summary="Normalize and annotate a single variant",
    description=(
        "Accepts one submitted variant in HGVS format, normalizes it through the "
        "ClinGen Allele Registry, annotates it through Ensembl VEP, and returns "
        "the internal AnnotatedVariant JSON shape. If annotation fails, the "
        "response still returns NormalizedVariant data plus annotation failure "
        "metadata."
    ),
)
def predict_single(
    variant: str = Query(
        description="Variant to normalize. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
) -> AnnotatedVariant:
    return _annotate_variant_or_raise(variant)


@router.get(
    "/annotateVariant",
    response_model=AnnotatedVariant,
    summary="Annotate a single variant",
    description=(
        "Accepts one submitted variant in HGVS format, normalizes it through the "
        "ClinGen Allele Registry, annotates it through Ensembl VEP, and returns "
        "the internal AnnotatedVariant JSON shape. If annotation fails, the "
        "response still returns NormalizedVariant data plus annotation failure "
        "metadata."
    ),
)
def annotate_single(
    variant: str = Query(
        description="Variant to annotate. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
) -> AnnotatedVariant:
    return _annotate_variant_or_raise(variant)


@router.post(
    "/predictOncogenicity",
    response_model=AnnotatedVariantBatchResponse,
    summary="Normalize and annotate a batch of variants",
    description=(
        "Accepts a list of submitted variants in HGVS format, normalizes each one "
        "through ClinGen, annotates each one through Ensembl VEP, and returns a "
        "list of internal AnnotatedVariant objects. Variants with annotation "
        "failures are returned in-band with annotation failure metadata rather "
        "than failing the whole batch."
    ),
)
def predict_batch(request: BatchRequest) -> AnnotatedVariantBatchResponse:
    annotated_variants = [
        _annotate_variant_or_raise(variant)
        for variant in request.variants
    ]
    return AnnotatedVariantBatchResponse(annotated_variants=annotated_variants)
