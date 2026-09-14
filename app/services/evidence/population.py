"""Population evidence scoring from summarized gnomAD frequencies.

This module consumes the collapsed population summary produced during
annotation rather than raw VEP frequency maps. It applies the current v1 rule
set for SBVS1, SBS1, and OP4 and preserves the key values that drove the
decision in matchedData.
"""

from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import EvidenceResult


def _format_percent(value: float) -> str:
    """Format a frequency as a percentage for clinician-facing evidence statements."""
    return f"{value * 100:.2f}%"


def build_population_evidence(annotated_variant: AnnotatedVariant) -> EvidenceResult:
    """Score the population pipeline from the annotation-layer frequency summary."""
    if annotated_variant.basicAnnotation is None:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "Population evidence could not be evaluated because annotation "
                "data was unavailable."
            ),
            status="not_available",
            source="vep",
            matchedData=None,
            dataAbsentReason="error",
        )

    population = annotated_variant.basicAnnotation.population
    effective_af_source = "none"
    effective_af = 0.0

    if population.maxSubpopulationAf is not None:
        effective_af = population.maxSubpopulationAf
        effective_af_source = "maxSubpopulationAf"
    elif population.maxOverallAf is not None:
        effective_af = population.maxOverallAf
        effective_af_source = "maxOverallAf"

    matched_data = {
        "effectiveAf": effective_af,
        "effectiveAfSource": effective_af_source,
        "maxSubpopulationAf": population.maxSubpopulationAf,
        "maxSubpopulationLabel": population.maxSubpopulationLabel,
        "maxOverallAf": population.maxOverallAf,
        "maxOverallLabel": population.maxOverallLabel,
    }

    if effective_af > 0.05:
        return EvidenceResult(
            score=-8,
            evidenceCode="SBVS1",
            evidenceStatement=(
                "Minor allele frequency is >5% in gnomAD "
                f"({_format_percent(effective_af)})."
            ),
            status="applied",
            source="vep",
            matchedData=matched_data,
        )

    if effective_af > 0.01:
        return EvidenceResult(
            score=-4,
            evidenceCode="SBS1",
            evidenceStatement=(
                "Minor allele frequency is >1% in gnomAD "
                f"({_format_percent(effective_af)})."
            ),
            status="applied",
            source="vep",
            matchedData=matched_data,
        )

    if effective_af == 0.0:
        statement = "Absent from gnomAD controls."
    else:
        statement = (
            "Present at low frequency in gnomAD (≤1%; "
            f"observed {_format_percent(effective_af)})."
        )

    return EvidenceResult(
        score=1,
        evidenceCode="OP4",
        evidenceStatement=statement,
        status="applied",
        source="vep",
        matchedData=matched_data,
    )
