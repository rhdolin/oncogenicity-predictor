from urllib.parse import quote

import httpx

from ...models.annotated_variant import (
    AnnotationError,
    AnnotatedVariant,
    BasicAnnotation,
    CaddAnnotation,
    ComputationalAnnotation,
    FathmmXfCodingAnnotation,
    PopulationSummary,
    TranscriptConsequence,
)
from ...models.normalized_variant import NormalizedVariant


VEP_HGVS_URL = "https://rest.ensembl.org/vep/human/hgvs"
SUBPOPULATION_KEYS = ("afr", "eas", "nfe", "amr", "sas")
INVALID_PREDICTOR_VALUE = "invalid_field"


class VariantAnnotationError(Exception):
    def __init__(self, message: str, attempted_queries: list[str]) -> None:
        super().__init__(message)
        self.message = message
        self.attempted_queries = attempted_queries


def annotate_variant(normalized_variant: NormalizedVariant) -> AnnotatedVariant:
    try:
        record = fetch_vep_annotation_record(normalized_variant)
    except VariantAnnotationError as exc:
        return AnnotatedVariant(
            normalizedVariant=normalized_variant,
            annotationStatus="failed",
            annotationError=AnnotationError(
                source="vep",
                message=exc.message,
                attemptedQueries=exc.attempted_queries,
            ),
            basicAnnotation=None,
            computationalAnnotation=None,
        )

    transcript_rows = record.get("transcript_consequences") or []
    _backfill_mane_from_vep(normalized_variant, transcript_rows)

    return AnnotatedVariant(
        normalizedVariant=normalized_variant,
        annotationStatus="complete",
        annotationError=None,
        basicAnnotation=BasicAnnotation(
            mostSevereConsequence=record.get("most_severe_consequence"),
            transcriptConsequences=_extract_transcript_consequences(transcript_rows),
            population=_extract_population_summary(record),
        ),
        computationalAnnotation=ComputationalAnnotation(
            cadd=_extract_cadd_annotation(transcript_rows),
            phyloP100wayVertebrate=_extract_first_predictor_value(
                transcript_rows,
                "phylop100way_vertebrate",
            ),
            fathmmXfCoding=_extract_fathmm_xf_coding_annotation(transcript_rows),
        ),
    )


def fetch_vep_annotation_record(normalized_variant: NormalizedVariant) -> dict:
    query_candidates = _build_vep_query_candidates(normalized_variant)

    for query_value in query_candidates:
        encoded_query = quote(query_value, safe="")
        try:
            response = httpx.get(
                f"{VEP_HGVS_URL}/{encoded_query}",
                params={
                    "content-type": "application/json",
                    "CADD": "true",
                    # Ensembl REST currently returns invalid_field for explicit
                    # hyphenated FATHMM-XF field requests, but dbNSFP=ALL yields
                    # the populated transcript keys we need.
                    "dbNSFP": "ALL",
                    "hgvs": "1",
                    "refseq": "true",
                    "mane": "1",
                },
                headers={"Accept": "application/json"},
                timeout=30.0,
            )
        except httpx.HTTPError:
            continue

        if response.status_code != 200:
            continue

        payload = response.json()
        if isinstance(payload, list) and payload and isinstance(payload[0], dict):
            return payload[0]

    raise VariantAnnotationError(
        "VEP annotation failed for all supported query forms.",
        query_candidates,
    )


def _build_vep_query_candidates(normalized_variant: NormalizedVariant) -> list[str]:
    candidates: list[str] = []
    if normalized_variant.genomic_hgvs.GRCh38:
        candidates.append(normalized_variant.genomic_hgvs.GRCh38)
    if normalized_variant.transcript_hgvs.mane_select_b38:
        candidates.append(normalized_variant.transcript_hgvs.mane_select_b38)
    return candidates


def _extract_transcript_consequences(transcript_rows: list[dict]) -> list[TranscriptConsequence]:
    transcript_consequences: list[TranscriptConsequence] = []

    for row in transcript_rows:
        transcript_refseq = row.get("transcript_id")
        if not isinstance(transcript_refseq, str) or not transcript_refseq.startswith("NM_"):
            continue

        consequence_terms = row.get("consequence_terms")
        if not isinstance(consequence_terms, list):
            consequence_terms = []

        transcript_consequences.append(
            TranscriptConsequence(
                transcriptRefSeq=transcript_refseq,
                consequenceTerms=[term for term in consequence_terms if isinstance(term, str)],
                proteinStart=row.get("protein_start"),
                proteinEnd=row.get("protein_end"),
                aminoAcids=row.get("amino_acids"),
                isManeSelect=_is_mane_select_row(row),
            )
        )

    return transcript_consequences


def _backfill_mane_from_vep(normalized_variant: NormalizedVariant, transcript_rows: list[dict]) -> None:
    if normalized_variant.transcript_hgvs.mane_select_b38:
        return

    for row in transcript_rows:
        if not _is_mane_select_row(row):
            continue

        mane_hgvsc = row.get("hgvsc")
        if isinstance(mane_hgvsc, str) and mane_hgvsc.startswith("NM_"):
            normalized_variant.transcript_hgvs.mane_select_b38 = mane_hgvsc
            normalized_variant.transcript_hgvs.mane_select_b38_source = "vep"
            return


def _is_mane_select_row(transcript_row: dict) -> bool | None:
    if transcript_row.get("mane_select"):
        return True
    mane = transcript_row.get("mane")
    if mane == "MANE Select":
        return True
    if isinstance(mane, list) and "MANE_Select" in mane:
        return True
    return None


def _extract_population_summary(record: dict) -> PopulationSummary:
    try:
        alt = (record.get("allele_string") or "").split("/")[1]
    except IndexError:
        alt = None

    max_subpopulation_af = None
    max_subpopulation_label = None
    max_overall_af = None
    max_overall_label = None

    for colocated_variant in record.get("colocated_variants") or []:
        frequencies = colocated_variant.get("frequencies") or {}
        if not frequencies:
            continue

        allele_key = alt if alt in frequencies else next(iter(frequencies), None)
        if not allele_key:
            continue

        allele_frequencies = frequencies.get(allele_key) or {}
        for source_prefix in ("gnomade", "gnomadg"):
            for population_key in SUBPOPULATION_KEYS:
                frequency_key = f"{source_prefix}_{population_key}"
                frequency_value = allele_frequencies.get(frequency_key)
                if isinstance(frequency_value, (int, float)) and (
                    max_subpopulation_af is None or frequency_value > max_subpopulation_af
                ):
                    max_subpopulation_af = float(frequency_value)
                    max_subpopulation_label = frequency_key

            overall_value = allele_frequencies.get(source_prefix)
            if isinstance(overall_value, (int, float)) and (
                max_overall_af is None or overall_value > max_overall_af
            ):
                max_overall_af = float(overall_value)
                max_overall_label = source_prefix

    return PopulationSummary(
        maxSubpopulationAf=max_subpopulation_af,
        maxSubpopulationLabel=max_subpopulation_label,
        maxOverallAf=max_overall_af,
        maxOverallLabel=max_overall_label,
    )


def _extract_cadd_annotation(transcript_rows: list[dict]) -> CaddAnnotation:
    return CaddAnnotation(
        phred=_extract_first_predictor_value(transcript_rows, "cadd_phred"),
        raw=_extract_first_predictor_value(transcript_rows, "cadd_raw"),
    )


def _extract_fathmm_xf_coding_annotation(transcript_rows: list[dict]) -> FathmmXfCodingAnnotation:
    return FathmmXfCodingAnnotation(
        prediction=_extract_first_predictor_string_value(transcript_rows, "fathmm-xf_coding_pred"),
        score=_extract_first_predictor_value(transcript_rows, "fathmm-xf_coding_score"),
        rankscore=_extract_first_predictor_value(transcript_rows, "fathmm-xf_coding_rankscore"),
    )


def _extract_first_predictor_value(transcript_rows: list[dict], key: str) -> float | None:
    for row in _iter_prioritized_transcript_rows(transcript_rows):
        value = row.get(key)
        if isinstance(value, (int, float)):
            return float(value)

    return None


def _extract_first_predictor_string_value(transcript_rows: list[dict], key: str) -> str | None:
    for row in _iter_prioritized_transcript_rows(transcript_rows):
        value = row.get(key)
        if isinstance(value, str) and value and value.lower() != INVALID_PREDICTOR_VALUE:
            return value

    return None


def _iter_prioritized_transcript_rows(transcript_rows: list[dict]):
    prioritized_rows = [
        [row for row in transcript_rows if _is_mane_select_row(row)],
        [
            row
            for row in transcript_rows
            if isinstance(row.get("transcript_id"), str) and row["transcript_id"].startswith("NM_")
        ],
        transcript_rows,
    ]

    for rows in prioritized_rows:
        for row in rows:
            yield row
