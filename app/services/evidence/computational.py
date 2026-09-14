from app.models.annotated_variant import AnnotatedVariant, ComputationalAnnotation
from app.models.prediction import EvidenceResult


OP1_CADD_PHRED_THRESHOLD = 15.0
SBP1_BENIGN_PREDICTIONS = {"N", "neutral", "benign", "tolerated"}


def build_computational_evidence(annotated_variant: AnnotatedVariant) -> EvidenceResult:
    if (
        annotated_variant.annotationStatus != "complete"
        or annotated_variant.basicAnnotation is None
        or annotated_variant.computationalAnnotation is None
    ):
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="Computational evidence could not be evaluated because annotation data was unavailable.",
            status="not_available",
            source="vep",
            matchedData=None,
            dataAbsentReason="error",
        )

    computational_annotation = annotated_variant.computationalAnnotation
    most_severe_consequence = annotated_variant.basicAnnotation.mostSevereConsequence

    if most_severe_consequence != "missense_variant":
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "Computational missense rules were not applicable because the most severe consequence "
                f"was {most_severe_consequence or 'unknown'}."
            ),
            status="applied",
            source="vep",
            matchedData=_build_matched_data(computational_annotation, most_severe_consequence),
            dataAbsentReason=None,
        )

    cadd_phred = computational_annotation.cadd.phred
    if cadd_phred is None:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="Computational evidence could not be evaluated because CADD annotation was unavailable.",
            status="not_available",
            source="vep",
            matchedData=_build_matched_data(computational_annotation, most_severe_consequence),
            dataAbsentReason="unknown",
        )

    if cadd_phred >= OP1_CADD_PHRED_THRESHOLD:
        return EvidenceResult(
            score=1,
            evidenceCode="OP1",
            evidenceStatement=f"CADD supports oncogenicity for this missense variant (PHRED {cadd_phred:.1f}).",
            status="applied",
            source="vep",
            matchedData=_build_matched_data(computational_annotation, most_severe_consequence),
            dataAbsentReason=None,
        )

    fathmm_prediction = computational_annotation.fathmmXfCoding.prediction
    if _is_benign_fathmm_xf_prediction(fathmm_prediction):
        return EvidenceResult(
            score=-1,
            evidenceCode="SBP1",
            evidenceStatement=(
                "Concordant computational predictors support a benign effect for this missense variant "
                f"(CADD PHRED {cadd_phred:.1f}; FATHMM-XF {fathmm_prediction})."
            ),
            status="applied",
            source="vep",
            matchedData=_build_matched_data(computational_annotation, most_severe_consequence),
            dataAbsentReason=None,
        )

    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement="Computational evidence did not meet current scoring criteria.",
        status="applied",
        source="vep",
        matchedData=_build_matched_data(computational_annotation, most_severe_consequence),
        dataAbsentReason=None,
    )


def _is_benign_fathmm_xf_prediction(prediction: str | None) -> bool:
    if prediction is None:
        return False
    return prediction.strip().lower() in {value.lower() for value in SBP1_BENIGN_PREDICTIONS}


def _build_matched_data(
    computational_annotation: ComputationalAnnotation,
    most_severe_consequence: str | None,
) -> dict[str, float | str | None]:
    return {
        "mostSevereConsequence": most_severe_consequence,
        "caddPhred": computational_annotation.cadd.phred,
        "caddRaw": computational_annotation.cadd.raw,
        "phyloP100wayVertebrate": computational_annotation.phyloP100wayVertebrate,
        "fathmmXfCodingPrediction": computational_annotation.fathmmXfCoding.prediction,
        "fathmmXfCodingScore": computational_annotation.fathmmXfCoding.score,
        "fathmmXfCodingRankscore": computational_annotation.fathmmXfCoding.rankscore,
    }
