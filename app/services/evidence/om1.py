"""OM1 evidence scoring for critical and well-established functional domains."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from app.models.annotated_variant import AnnotatedVariant, TranscriptConsequence
from app.models.prediction import EvidenceResult


REPO_ROOT = Path(__file__).resolve().parents[3]
OM1_DOMAINS_PATH = REPO_ROOT / "data" / "om1_clingen_domains_seed.csv"
ELIGIBLE_EVENT_TYPES = {"substitution", "deletion", "insertion", "duplication", "delins"}


class Om1DomainsError(Exception):
    """Raised when the local OM1 domains table cannot be parsed."""


@dataclass(frozen=True, slots=True)
class Om1DomainRule:
    gene_symbol: str
    mane_transcript: str
    mane_protein: str | None
    domain_name: str
    start_residue: int
    end_residue: int
    excluded_residues: frozenset[int]
    source: str | None
    source_version: str | None
    notes: str | None


def build_om1_evidence(annotated_variant: AnnotatedVariant) -> EvidenceResult:
    """Evaluate OM1 from the curated ClinGen domain table."""
    if annotated_variant.annotationStatus != "complete" or annotated_variant.basicAnnotation is None:
        return _build_not_available_result(
            statement="OM1 evidence could not be evaluated because annotation data was unavailable.",
            matched_data=None,
            absent_reason="error",
        )

    gene_symbol = (annotated_variant.normalizedVariant.geneSymbol or "").upper()
    matched_data = {"gene": gene_symbol or None}

    if not gene_symbol:
        return _build_not_available_result(
            statement="OM1 evidence could not be evaluated because the normalized gene symbol was unavailable.",
            matched_data=matched_data,
            absent_reason="unknown",
        )

    try:
        gene_rows, ready_rows = _get_gene_domain_rows(gene_symbol)
    except (FileNotFoundError, Om1DomainsError) as exc:
        return _build_not_available_result(
            statement=(
                "OM1 evidence could not be evaluated because the local OM1 domains "
                f"table was unavailable or unreadable: {exc}"
            ),
            matched_data=matched_data,
            absent_reason="error",
        )

    matched_data["geneRowCount"] = len(gene_rows)
    matched_data["readyDomainCount"] = len(ready_rows)

    if not gene_rows:
        return _build_not_available_result(
            statement=(
                "OM1 evidence not available: gene "
                f"{gene_symbol} is not represented in the local ClinGen domain table."
            ),
            matched_data=matched_data,
        )

    if not ready_rows:
        return _build_not_available_result(
            statement=(
                "OM1 evidence not available: gene "
                f"{gene_symbol} has no released ClinGen critical-domain rows in the local table."
            ),
            matched_data=matched_data,
        )

    mane_consequence = _extract_mane_consequence(
        annotated_variant.basicAnnotation.transcriptConsequences
    )
    if mane_consequence is None:
        return _build_not_available_result(
            statement=(
                "OM1 evidence not available: a MANE Select transcript consequence "
                f"was required for gene {gene_symbol}."
            ),
            matched_data=matched_data,
        )

    matched_data.update(
        {
            "maneTranscript": mane_consequence.transcriptRefSeq,
            "proteinHgvs": mane_consequence.proteinHgvs,
            "proteinEventType": mane_consequence.proteinEventType,
            "proteinStart": mane_consequence.proteinStart,
            "proteinEnd": mane_consequence.proteinEnd,
        }
    )

    if mane_consequence.proteinEventType not in ELIGIBLE_EVENT_TYPES:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="OM1 evidence did not meet current scoring criteria.",
            status="applied",
            source="clingen-cspec",
            matchedData=matched_data,
            dataAbsentReason=None,
        )

    compatible_rows = [
        row for row in ready_rows if row.mane_transcript == mane_consequence.transcriptRefSeq
    ]
    matched_data["compatibleDomainCount"] = len(compatible_rows)

    if not compatible_rows:
        return _build_not_available_result(
            statement=(
                "OM1 evidence not available: no released ClinGen domain definitions "
                f"matched MANE transcript {mane_consequence.transcriptRefSeq} for gene {gene_symbol}."
            ),
            matched_data=matched_data,
        )

    event_start, event_end = _resolve_event_bounds(mane_consequence)
    if event_start is None or event_end is None:
        return _build_not_available_result(
            statement=(
                "OM1 evidence not available: a localized MANE protein residue position "
                f"or span could not be resolved for gene {gene_symbol}."
            ),
            matched_data=matched_data,
        )

    matched_data["eventStart"] = event_start
    matched_data["eventEnd"] = event_end

    for row in compatible_rows:
        if not _intervals_overlap(event_start, event_end, row.start_residue, row.end_residue):
            continue
        excluded_overlap = sorted(
            residue
            for residue in row.excluded_residues
            if event_start <= residue <= event_end
        )
        if excluded_overlap:
            continue

        return EvidenceResult(
            score=2,
            evidenceCode="OM1",
            evidenceStatement=(
                "Located in a critical and well-established functional domain "
                "defined in the curated local ClinGen domain table."
            ),
            status="applied",
            source="clingen-cspec",
            matchedData={
                **matched_data,
                "matchedDomain": {
                    "geneSymbol": row.gene_symbol,
                    "maneTranscript": row.mane_transcript,
                    "maneProtein": row.mane_protein,
                    "domainName": row.domain_name,
                    "startResidue": row.start_residue,
                    "endResidue": row.end_residue,
                    "excludedResidues": sorted(row.excluded_residues),
                    "source": row.source,
                    "sourceVersion": row.source_version,
                    "notes": row.notes,
                },
            },
            dataAbsentReason=None,
        )

    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement="OM1 evidence did not meet current scoring criteria.",
        status="applied",
        source="clingen-cspec",
        matchedData=matched_data,
        dataAbsentReason=None,
    )


def _build_not_available_result(
    statement: str,
    matched_data: dict | None,
    absent_reason: str = "unsupported",
) -> EvidenceResult:
    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement=statement,
        status="not_available",
        source="clingen-cspec",
        matchedData=matched_data,
        dataAbsentReason=absent_reason,
    )


def _extract_mane_consequence(
    transcript_consequences: list[TranscriptConsequence],
) -> TranscriptConsequence | None:
    for consequence in transcript_consequences:
        if consequence.isManeSelect:
            return consequence
    return None


def _resolve_event_bounds(
    consequence: TranscriptConsequence,
) -> tuple[int | None, int | None]:
    if consequence.proteinEventType == "substitution":
        if consequence.proteinStart is not None and consequence.proteinEnd is not None:
            return consequence.proteinStart, consequence.proteinEnd
        if not consequence.proteinHgvs:
            return None, None
        match = re.match(r"^p\.[A-Z*](\d+)[A-Z*=]$", consequence.proteinHgvs)
        if not match:
            return None, None
        position = int(match.group(1))
        return position, position

    if consequence.proteinHgvs:
        event_start, event_end = _parse_localized_protein_event_bounds(
            consequence.proteinHgvs.removeprefix("p.")
        )
        if event_start is not None and event_end is not None:
            return event_start, event_end

    return consequence.proteinStart, consequence.proteinEnd


def _intervals_overlap(
    left_start: int,
    left_end: int,
    right_start: int,
    right_end: int,
) -> bool:
    return left_start <= right_end and right_start <= left_end


def _parse_localized_protein_event_bounds(event: str) -> tuple[int | None, int | None]:
    # Hotspots uses the same local event-shape parser.
    # If these rules change, review app/services/evidence/hotspots.py as well.
    match = re.match(
        r"^[A-Z*](\d+)(?:_[A-Z*](\d+))?(?:delins[A-Z*]+|ins[A-Z*]+|del|dup)$",
        event,
    )
    if not match:
        return None, None

    start_text, end_text = match.groups()
    start = int(start_text)
    end = int(end_text) if end_text is not None else start
    return start, end


def _parse_excluded_residues(value: str | None) -> frozenset[int]:
    if value is None:
        return frozenset()

    residues: set[int] = set()
    for token in value.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        range_match = re.match(r"^(\d+)-(\d+)$", stripped)
        if range_match:
            start_text, end_text = range_match.groups()
            start = int(start_text)
            end = int(end_text)
            residues.update(range(start, end + 1))
            continue
        if stripped.isdigit():
            residues.add(int(stripped))
    return frozenset(residues)


@lru_cache(maxsize=1)
def _load_om1_domain_rows() -> tuple[Om1DomainRule, ...]:
    if not OM1_DOMAINS_PATH.exists():
        raise FileNotFoundError(OM1_DOMAINS_PATH)

    rows: list[Om1DomainRule] = []
    with OM1_DOMAINS_PATH.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        expected_columns = {
            "geneSymbol",
            "maneTranscript",
            "maneProtein",
            "domainName",
            "startResidue",
            "endResidue",
            "excludedResidues",
            "rowStatus",
            "source",
            "sourceVersion",
            "notes",
        }
        if set(reader.fieldnames or []) != expected_columns:
            raise Om1DomainsError(
                "expected columns geneSymbol,maneTranscript,maneProtein,"
                "domainName,startResidue,endResidue,excludedResidues,rowStatus,"
                "source,sourceVersion,notes"
            )

        for row in reader:
            if (row.get("rowStatus") or "").strip() != "ready":
                continue

            gene_symbol = (row.get("geneSymbol") or "").strip().upper()
            mane_transcript = (row.get("maneTranscript") or "").strip()
            domain_name = (row.get("domainName") or "").strip()
            start_residue = _parse_required_int(row.get("startResidue"))
            end_residue = _parse_required_int(row.get("endResidue"))
            if not gene_symbol or not mane_transcript or not domain_name:
                raise Om1DomainsError(
                    "ready OM1 rows must include geneSymbol, maneTranscript, "
                    "and domainName"
                )
            rows.append(
                Om1DomainRule(
                    gene_symbol=gene_symbol,
                    mane_transcript=mane_transcript,
                    mane_protein=(row.get("maneProtein") or "").strip() or None,
                    domain_name=domain_name,
                    start_residue=start_residue,
                    end_residue=end_residue,
                    excluded_residues=_parse_excluded_residues(row.get("excludedResidues")),
                    source=(row.get("source") or "").strip() or None,
                    source_version=(row.get("sourceVersion") or "").strip() or None,
                    notes=(row.get("notes") or "").strip() or None,
                )
            )

    return tuple(rows)


@lru_cache(maxsize=1)
def _get_gene_domain_rows(
    gene_symbol: str,
) -> tuple[tuple[dict[str, str], ...], tuple[Om1DomainRule, ...]]:
    ready_rows = tuple(row for row in _load_om1_domain_rows() if row.gene_symbol == gene_symbol)

    if not OM1_DOMAINS_PATH.exists():
        raise FileNotFoundError(OM1_DOMAINS_PATH)

    all_gene_rows: list[dict[str, str]] = []
    with OM1_DOMAINS_PATH.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        expected_columns = {
            "geneSymbol",
            "maneTranscript",
            "maneProtein",
            "domainName",
            "startResidue",
            "endResidue",
            "excludedResidues",
            "rowStatus",
            "source",
            "sourceVersion",
            "notes",
        }
        if set(reader.fieldnames or []) != expected_columns:
            raise Om1DomainsError(
                "expected columns geneSymbol,maneTranscript,maneProtein,"
                "domainName,startResidue,endResidue,excludedResidues,rowStatus,"
                "source,sourceVersion,notes"
            )
        for row in reader:
            if (row.get("geneSymbol") or "").strip().upper() == gene_symbol:
                all_gene_rows.append({key: value or "" for key, value in row.items()})

    return tuple(all_gene_rows), ready_rows


def _parse_required_int(value: str | None) -> int:
    if value is None:
        raise Om1DomainsError("ready OM1 rows must include numeric residue bounds")
    stripped = value.strip()
    if not stripped:
        raise Om1DomainsError("ready OM1 rows must include numeric residue bounds")
    try:
        return int(stripped)
    except ValueError as exc:
        raise Om1DomainsError("ready OM1 rows must include numeric residue bounds") from exc
