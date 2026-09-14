"""HTTP route definitions for annotation, evidence summary, and prediction.

This module keeps transport concerns thin: request parsing, response models,
and translation of normalization failures into HTTP errors. Domain workflow is
deferred to orchestration services so the pipeline logic stays testable.
"""

from fastapi import APIRouter, HTTPException, Query

from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import (
    OncogenicityObservation,
    OncogenicityPredictionBatchResponse,
    OncogenicityPredictionSummary,
)
from app.models.requests import BatchRequest
from app.services.normalization.variant_normalizer import VariantNormalizationError
from app.services.orchestration.single_variant_pipeline import (
    run_single_variant_annotation_pipeline,
    run_single_variant_evidence_summary_pipeline,
    run_single_variant_pipeline,
)


router = APIRouter()


def _annotate_variant_or_raise(submitted_variant: str) -> AnnotatedVariant:
    """Run the annotation pipeline and translate normalization failures into HTTP 400 errors."""
    try:
        return run_single_variant_annotation_pipeline(submitted_variant)
    except VariantNormalizationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _predict_variant_or_raise(submitted_variant: str) -> OncogenicityObservation:
    """Run the FHIR prediction pipeline and translate normalization failures into HTTP 400 errors."""
    try:
        return run_single_variant_pipeline(submitted_variant)
    except VariantNormalizationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _summarize_variant_or_raise(submitted_variant: str) -> OncogenicityPredictionSummary:
    """Run the internal evidence-summary pipeline and translate normalization failures into HTTP 400 errors."""
    try:
        return run_single_variant_evidence_summary_pipeline(submitted_variant)
    except VariantNormalizationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get(
    "/predictOncogenicity",
    response_model=OncogenicityObservation,
    response_model_exclude_none=True,
    summary="Predict oncogenicity for a single variant",
    description=(
        "Accepts one submitted variant in HGVS format, normalizes it through the "
        "ClinGen Allele Registry, annotates it through Ensembl VEP, evaluates "
        "the currently implemented evidence pipelines, and returns a single "
        "FHIR Observation-style prediction object with `issued`, a custom "
        "extension carrying the submitted variant, an overall score, "
        "and one component per evidence pipeline."
    ),
)
def predict_single(
    variant: str = Query(
        description="Variant to normalize. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
) -> OncogenicityObservation:
    """Handle GET requests for a single clinician-facing FHIR prediction."""
    return _predict_variant_or_raise(variant)


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
    """Handle GET requests for the internal annotation payload."""
    return _annotate_variant_or_raise(variant)


@router.get(
    "/summarizeEvidence",
    response_model=OncogenicityPredictionSummary,
    summary="Summarize evidence for a single variant",
    description=(
        "Internal/debug endpoint. Accepts one submitted variant in HGVS format, "
        "normalizes it through the ClinGen Allele Registry, annotates it through "
        "Ensembl VEP, evaluates the currently implemented evidence pipelines, "
        "and returns the raw evidence summary JSON before FHIR Observation mapping."
    ),
)
def summarize_single(
    variant: str = Query(
        description="Variant to summarize. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
) -> OncogenicityPredictionSummary:
    """Handle GET requests for the raw evidence summary before FHIR mapping."""
    return _summarize_variant_or_raise(variant)


@router.post(
    "/predictOncogenicity",
    response_model=OncogenicityPredictionBatchResponse,
    response_model_exclude_none=True,
    summary="Predict oncogenicity for a batch of variants",
    description=(
        "Accepts a list of submitted variants in HGVS format, normalizes each one "
        "through ClinGen, annotates each one through Ensembl VEP, evaluates the "
        "currently implemented evidence pipelines, and returns a list of "
        "FHIR Observation-style prediction objects. Variants with annotation "
        "failures are returned in-band as partial observations with pipeline "
        "`dataAbsentReason` values rather than failing the whole batch."
    ),
)
def predict_batch(request: BatchRequest) -> OncogenicityPredictionBatchResponse:
    """Handle batch prediction requests by running the single-variant pipeline per input."""
    observations = [_predict_variant_or_raise(variant) for variant in request.variants]
    return OncogenicityPredictionBatchResponse(observations=observations)
