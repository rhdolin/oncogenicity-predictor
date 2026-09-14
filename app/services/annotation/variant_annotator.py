"""Annotation service backed by Ensembl REST VEP.

This module converts a normalized variant into the internal AnnotatedVariant
shape used by downstream evidence builders. It also owns VEP query selection,
lightweight transcript prioritization, and extraction of lean population and
computational annotations for scoring.
"""

import re
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
AMINO_ACID_THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}


class VariantAnnotationError(Exception):
    """Raised when VEP annotation fails for every supported query form."""

    def __init__(self, message: str, attempted_queries: list[str]) -> None:
        super().__init__(message)
        self.message = message
        self.attempted_queries = attempted_queries


def annotate_variant(normalized_variant: NormalizedVariant) -> AnnotatedVariant:
    """Annotate one normalized variant and return a partial failure payload if VEP cannot resolve it."""
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
    """Query VEP using the supported HGVS fallbacks and return the first usable record."""
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
    """Return the ordered list of HGVS query forms the current VEP strategy will try."""
    candidates: list[str] = []
    if normalized_variant.genomic_hgvs.GRCh38:
        candidates.append(normalized_variant.genomic_hgvs.GRCh38)
    if normalized_variant.transcript_hgvs.mane_select_b38:
        candidates.append(normalized_variant.transcript_hgvs.mane_select_b38)
    return candidates


def _extract_transcript_consequences(transcript_rows: list[dict]) -> list[TranscriptConsequence]:
    """Keep only RefSeq transcript rows and map them into the lean internal transcript model."""
    transcript_consequences: list[TranscriptConsequence] = []

    for row in transcript_rows:
        transcript_refseq = row.get("transcript_id")
        if not isinstance(transcript_refseq, str) or not transcript_refseq.startswith("NM_"):
            continue

        consequence_terms = row.get("consequence_terms")
        if not isinstance(consequence_terms, list):
            consequence_terms = []

        protein_hgvs = _extract_protein_hgvs(row)

        transcript_consequences.append(
            TranscriptConsequence(
                transcriptRefSeq=transcript_refseq,
                consequenceTerms=[term for term in consequence_terms if isinstance(term, str)],
                proteinStart=row.get("protein_start"),
                proteinEnd=row.get("protein_end"),
                aminoAcids=row.get("amino_acids"),
                proteinHgvs=protein_hgvs,
                proteinEventType=_extract_protein_event_type(row, protein_hgvs),
                rawProteinHgvs=row.get("hgvsp") if isinstance(row.get("hgvsp"), str) else None,
                isManeSelect=_is_mane_select_row(row),
            )
        )

    return transcript_consequences


def _backfill_mane_from_vep(normalized_variant: NormalizedVariant, transcript_rows: list[dict]) -> None:
    """Populate MANE Select transcript HGVS from VEP when ClinGen did not provide it."""
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


def _extract_protein_hgvs(transcript_row: dict) -> str | None:
    """Return the transcript-row protein HGVS in stripped one-letter form when VEP provides it."""
    protein_hgvs = transcript_row.get("hgvsp")
    if not isinstance(protein_hgvs, str) or not protein_hgvs:
        return None

    stripped_hgvs = protein_hgvs.split(":", maxsplit=1)[-1]
    return _convert_protein_hgvs_to_one_letter(stripped_hgvs)


def _extract_protein_event_type(transcript_row: dict, protein_hgvs: str | None) -> str | None:
    """Classify the transcript-row protein effect into a compact event type."""
    raw_protein_hgvs = transcript_row.get("hgvsp")

    if protein_hgvs:
        substitution_match = re.match(r"^p\.([A-Z*])(\d+)([A-Z*=])$", protein_hgvs)
        if substitution_match:
            ref_residue, _start_residue, alt_residue = substitution_match.groups()
            event_type = "substitution"
            if alt_residue == "=":
                event_type = "silent"
            elif alt_residue == "*":
                event_type = "stop_gain"
            elif ref_residue == "*":
                event_type = "stop_loss"
            return event_type

        if re.match(r"^p\.([A-Z*])(\d+)(?:_([A-Z*])(\d+))?delins([A-Z*]+)$", protein_hgvs):
            return "delins"

        if re.match(r"^p\.([A-Z*])(\d+)_([A-Z*])(\d+)ins([A-Z*]+)$", protein_hgvs):
            return "insertion"

        if re.match(r"^p\.([A-Z*])(\d+)(?:_([A-Z*])(\d+))?del$", protein_hgvs):
            return "deletion"

        if re.match(r"^p\.([A-Z*])(\d+)(?:_([A-Z*])(\d+))?dup$", protein_hgvs):
            return "duplication"

        frameshift_match = re.match(r"^p\.([A-Z*])(\d+)(?:[A-Z*])?fs", protein_hgvs)
        if frameshift_match:
            return "frameshift"

        if "ext" in protein_hgvs:
            return "extension"

    if not isinstance(raw_protein_hgvs, str) or not raw_protein_hgvs:
        return None

    return "unknown"


def _convert_protein_hgvs_to_one_letter(hgvs_value: str) -> str:
    def replace_match(match: re.Match[str]) -> str:
        return AMINO_ACID_THREE_TO_ONE.get(match.group(0), match.group(0))

    return re.sub(
        r"Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter",
        replace_match,
        hgvs_value,
    )


def _is_mane_select_row(transcript_row: dict) -> bool | None:
    """Detect whether a transcript row carries a MANE Select marker in the observed VEP shapes."""
    if transcript_row.get("mane_select"):
        return True
    mane = transcript_row.get("mane")
    if mane == "MANE Select":
        return True
    if isinstance(mane, list) and "MANE_Select" in mane:
        return True
    return None


def _extract_population_summary(record: dict) -> PopulationSummary:
    """Collapse raw VEP colocated-variant frequencies into the summary fields used by scoring."""
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
    """Extract the prioritized CADD values from transcript consequence rows."""
    return CaddAnnotation(
        phred=_extract_first_predictor_value(transcript_rows, "cadd_phred"),
        raw=_extract_first_predictor_value(transcript_rows, "cadd_raw"),
    )


def _extract_fathmm_xf_coding_annotation(transcript_rows: list[dict]) -> FathmmXfCodingAnnotation:
    """Extract the prioritized FATHMM-XF fields used by the computational evidence rule."""
    return FathmmXfCodingAnnotation(
        prediction=_extract_first_predictor_string_value(transcript_rows, "fathmm-xf_coding_pred"),
        score=_extract_first_predictor_value(transcript_rows, "fathmm-xf_coding_score"),
        rankscore=_extract_first_predictor_value(transcript_rows, "fathmm-xf_coding_rankscore"),
    )


def _extract_first_predictor_value(transcript_rows: list[dict], key: str) -> float | None:
    """Return the first numeric predictor value from the prioritized transcript row order."""
    for row in _iter_prioritized_transcript_rows(transcript_rows):
        value = row.get(key)
        if isinstance(value, (int, float)):
            return float(value)

    return None


def _extract_first_predictor_string_value(transcript_rows: list[dict], key: str) -> str | None:
    """Return the first non-empty predictor string value, ignoring known VEP sentinel values."""
    for row in _iter_prioritized_transcript_rows(transcript_rows):
        value = row.get(key)
        if isinstance(value, str) and value and value.lower() != INVALID_PREDICTOR_VALUE:
            return value

    return None


def _iter_prioritized_transcript_rows(transcript_rows: list[dict]):
    """Yield transcript rows in the precedence order used for computational field extraction."""
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
