"""HTTP route definitions for annotation, evidence summary, and prediction.

This module keeps transport concerns thin: request parsing, response models,
and translation of normalization failures into HTTP errors. Domain workflow is
deferred to orchestration services so the pipeline logic stays testable.
"""

from fastapi import APIRouter, HTTPException, Query

from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import (
    OncogenicityObservation,
    OncogenicityObservationBundle,
    OncogenicityPredictionSummary,
)
from app.models.requests import BatchRequest
from app.models.tumor_types import TUMOR_TYPE_VALUES, TumorType
from app.services.fhir import build_oncogenicity_observation_bundle
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


def _predict_variant_or_raise(
    submitted_variant: str,
    tumor_type: TumorType | None = None,
) -> OncogenicityObservation:
    """Run the FHIR prediction pipeline and translate normalization failures into HTTP 400 errors."""
    try:
        return run_single_variant_pipeline(submitted_variant, tumor_type=tumor_type)
    except VariantNormalizationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _summarize_variant_or_raise(
    submitted_variant: str,
    tumor_type: TumorType | None = None,
) -> OncogenicityPredictionSummary:
    """Run the internal evidence-summary pipeline and translate normalization failures into HTTP 400 errors."""
    try:
        return run_single_variant_evidence_summary_pipeline(
            submitted_variant,
            tumor_type=tumor_type,
        )
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
        "extension carrying the submitted variant, an overall score, a final "
        "classification when available, and score-contributing evidence components. "
        "If no evidence lanes are evaluable, the Observation is returned with "
        "a top-level data-absent reason instead. Also accepts an optional "
        "`tumorType` query parameter used by context-dependent evidence rules."
    ),
)
def predict_single(
    variant: str = Query(
        description="Variant for which to predict oncogenicity. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
    tumorType: TumorType | None = Query(
        default=None,
        description="Optional tumor type used for context-dependent evidence rules.",
        examples=list(TUMOR_TYPE_VALUES),
    ),
) -> OncogenicityObservation:
    """Handle GET requests for a single clinician-facing FHIR prediction."""
    return _predict_variant_or_raise(variant, tumor_type=tumorType)


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
        "and returns the raw evidence summary JSON before FHIR Observation mapping. "
        "If no evidence lanes are evaluable, the summary exposes a top-level "
        "data-absent reason and no final score or classification. "
        "Also accepts an optional `tumorType` query parameter used by context-dependent evidence rules."
    ),
)
def summarize_single(
    variant: str = Query(
        description="Variant to summarize. Must be provided in HGVS format.",
        examples=["NM_004119.3:c.2073T>G"],
    ),
    tumorType: TumorType | None = Query(
        default=None,
        description="Optional tumor type used for context-dependent evidence rules.",
        examples=list(TUMOR_TYPE_VALUES),
    ),
) -> OncogenicityPredictionSummary:
    """Handle GET requests for the raw evidence summary before FHIR mapping."""
    return _summarize_variant_or_raise(variant, tumor_type=tumorType)


@router.post(
    "/predictOncogenicity",
    response_model=OncogenicityObservationBundle,
    response_model_exclude_none=True,
    summary="Predict oncogenicity for a batch of variants",
    description=(
        "Accepts a list of submitted variants in HGVS format, normalizes each one "
        "through ClinGen, annotates each one through Ensembl VEP, evaluates the "
        "currently implemented evidence pipelines, and returns a FHIR Bundle "
        "containing one Observation-style prediction resource per input variant. Variants with annotation "
        "failures are returned in-band rather than failing the whole batch, but "
        "the clinician-facing FHIR component list remains compact and may "
        "therefore be empty for those variants, with a top-level data-absent "
        "reason when the overall prediction is unavailable. The request body "
        "also accepts an optional top-level `tumorType` field applied to each batch entry."
    ),
)
def predict_batch(request: BatchRequest) -> OncogenicityObservationBundle:
    """Handle batch prediction requests by running the single-variant pipeline per input."""
    observations = [
        _predict_variant_or_raise(variant, tumor_type=request.tumorType)
        for variant in request.variants
    ]
    return build_oncogenicity_observation_bundle(observations)
