"""Orchestration entrypoints for the single-variant workflow.

These functions connect the major pipeline stages without embedding HTTP or
source-specific transport logic. They define the sequence from normalization to
annotation to either raw evidence summary generation or FHIR projection.
"""

from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import OncogenicityObservation, OncogenicityPredictionSummary
from app.models.tumor_types import TumorType
from app.services.annotation.variant_annotator import annotate_variant
from app.services.fhir import build_oncogenicity_observation
from app.services.normalization.variant_normalizer import normalize_variant

from .prediction_builder import build_prediction_summary_from_annotated_variant


def run_single_variant_annotation_pipeline(submitted_variant: str) -> AnnotatedVariant:
    """Normalize one submitted variant and return its annotation-layer payload."""
    normalized_variant = normalize_variant(submitted_variant)
    return annotate_variant(normalized_variant)


def run_single_variant_pipeline(
    submitted_variant: str,
    tumor_type: TumorType | None = None,
) -> OncogenicityObservation:
    """Run the full single-variant workflow and project the result into FHIR Observation form."""
    summary = run_single_variant_evidence_summary_pipeline(
        submitted_variant,
        tumor_type=tumor_type,
    )
    return build_oncogenicity_observation(summary, submitted_variant)


def run_single_variant_evidence_summary_pipeline(
    submitted_variant: str,
    tumor_type: TumorType | None = None,
) -> OncogenicityPredictionSummary:
    """Run the single-variant workflow through evidence scoring but stop before FHIR mapping."""
    annotated_variant = run_single_variant_annotation_pipeline(submitted_variant)
    summary = build_prediction_summary_from_annotated_variant(
        annotated_variant,
        tumor_type=tumor_type,
    )
    return summary
