"""Cancer Hotspots evidence scoring.

This module loads the bundled Cancer Hotspots workbook into an in-memory index
on first use and applies the current hotspot rule set against transcript-level
protein consequences. SNVs require exact gene, position, and amino-acid change
agreement. Indels are supported only when the transcript consequence can be
compared directly to the workbook's native event representation.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re
from zipfile import ZipFile
import xml.etree.ElementTree as ET

from app.models.annotated_variant import AnnotatedVariant, TranscriptConsequence
from app.models.prediction import EvidenceResult


WORKBOOK_PATH = Path(__file__).resolve().parents[3] / "data" / "hotspots_v2.xlsx"
XML_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
INDIRECT_EVENT_TYPES = {"deletion", "insertion", "duplication", "delins"}


class HotspotDataError(Exception):
    """Raised when the local hotspot workbook cannot be parsed into a usable index."""


@dataclass(frozen=True)
class SnpHotspotRecord:
    gene: str
    position: int
    ref_amino_acid: str
    alt_amino_acid: str
    mutation_count: int
    variant_count: int


@dataclass(frozen=True)
class IndelHotspotRecord:
    gene: str
    hotspot_position: str
    event: str
    event_start: int | None
    event_end: int | None
    mutation_count: int
    variant_count: int


@dataclass(frozen=True)
class HotspotIndex:
    snv_by_gene: dict[str, tuple[SnpHotspotRecord, ...]]
    indel_by_gene: dict[str, tuple[IndelHotspotRecord, ...]]


def build_hotspots_evidence(annotated_variant: AnnotatedVariant) -> EvidenceResult:
    """Score Cancer Hotspots evidence for one annotated variant."""
    if annotated_variant.annotationStatus != "complete" or annotated_variant.basicAnnotation is None:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="Hotspots evidence could not be evaluated because annotation data was unavailable.",
            status="not_available",
            source="cancerhotspots",
            matchedData=None,
            dataAbsentReason="error",
        )

    gene_symbol = annotated_variant.normalizedVariant.geneSymbol
    transcript_consequences = annotated_variant.basicAnnotation.transcriptConsequences

    if not gene_symbol:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "Hotspots evidence could not be evaluated because the normalized "
                "gene symbol was unavailable."
            ),
            status="not_available",
            source="cancerhotspots",
            matchedData=None,
            dataAbsentReason="unknown",
        )

    try:
        hotspot_index = _get_hotspot_index()
    except (FileNotFoundError, HotspotDataError) as exc:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "Hotspots evidence could not be evaluated because the local Cancer Hotspots data "
                f"was unavailable or unreadable: {exc}"
            ),
            status="not_available",
            source="cancerhotspots",
            matchedData={"gene": gene_symbol},
            dataAbsentReason="error",
        )

    for consequence in transcript_consequences:
        matched_record = _match_transcript_consequence(gene_symbol, consequence, hotspot_index)
        if matched_record is not None:
            return matched_record

    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement="Hotspots evidence did not meet current scoring criteria.",
        status="applied",
        source="cancerhotspots",
        matchedData={
            "gene": gene_symbol,
            "evaluatedTranscriptCount": len(transcript_consequences),
        },
        dataAbsentReason=None,
    )


@lru_cache(maxsize=1)
def _get_hotspot_index() -> HotspotIndex:
    workbook_rows = _read_workbook_rows(WORKBOOK_PATH)
    return HotspotIndex(
        snv_by_gene=_build_snv_index(workbook_rows.get("SNV-hotspots", [])),
        indel_by_gene=_build_indel_index(workbook_rows.get("INDEL-hotspots", [])),
    )


def _read_workbook_rows(workbook_path: Path) -> dict[str, list[dict[str, str]]]:
    if not workbook_path.exists():
        raise FileNotFoundError(workbook_path)

    with ZipFile(workbook_path) as workbook_zip:
        shared_strings = _load_shared_strings(workbook_zip)
        sheet_targets = _load_sheet_targets(workbook_zip)
        rows_by_sheet: dict[str, list[dict[str, str]]] = {}

        for sheet_name, target in sheet_targets.items():
            rows = list(_iter_sheet_rows(workbook_zip, target, shared_strings))
            if not rows:
                rows_by_sheet[sheet_name] = []
                continue

            header_row = rows[0]
            headers_by_column = {
                column: value.strip()
                for column, value in header_row.items()
                if value and value.strip()
            }

            data_rows: list[dict[str, str]] = []
            for row in rows[1:]:
                mapped_row = {
                    header: row.get(column, "").strip()
                    for column, header in headers_by_column.items()
                }
                if any(mapped_row.values()):
                    data_rows.append(mapped_row)

            rows_by_sheet[sheet_name] = data_rows

    return rows_by_sheet


def _load_shared_strings(workbook_zip: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in workbook_zip.namelist():
        return []

    root = ET.fromstring(workbook_zip.read("xl/sharedStrings.xml"))
    shared_strings: list[str] = []

    for item in root.findall("a:si", XML_NS):
        text_parts = [text_node.text or "" for text_node in item.iterfind(".//a:t", XML_NS)]
        shared_strings.append("".join(text_parts))

    return shared_strings


def _load_sheet_targets(workbook_zip: ZipFile) -> dict[str, str]:
    workbook_root = ET.fromstring(workbook_zip.read("xl/workbook.xml"))
    relationships_root = ET.fromstring(workbook_zip.read("xl/_rels/workbook.xml.rels"))

    relationships = {
        relationship.attrib.get("Id"): relationship.attrib.get("Target", "")
        for relationship in relationships_root
    }

    sheet_targets: dict[str, str] = {}
    for sheet in workbook_root.findall("a:sheets/a:sheet", XML_NS):
        relationship_id = sheet.attrib.get(
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        )
        sheet_name = sheet.attrib.get("name")
        if not relationship_id or not sheet_name:
            continue
        target = relationships.get(relationship_id)
        if target:
            sheet_targets[sheet_name] = target

    return sheet_targets


def _iter_sheet_rows(
    workbook_zip: ZipFile,
    sheet_target: str,
    shared_strings: list[str],
) -> list[dict[str, str]]:
    sheet_root = ET.fromstring(workbook_zip.read("xl/" + sheet_target.lstrip("/")))

    for row in sheet_root.findall("a:sheetData/a:row", XML_NS):
        row_values: dict[str, str] = {}
        for cell in row.findall("a:c", XML_NS):
            cell_reference = cell.attrib.get("r", "")
            column_name = "".join(character for character in cell_reference if character.isalpha())
            row_values[column_name] = _extract_cell_text(cell, shared_strings)
        yield row_values


def _extract_cell_text(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(text_node.text or "" for text_node in cell.iterfind(".//a:t", XML_NS))

    raw_value = cell.findtext("a:v", default="", namespaces=XML_NS)
    if cell_type == "s" and raw_value:
        return shared_strings[int(raw_value)]
    return raw_value


def _build_snv_index(rows: list[dict[str, str]]) -> dict[str, tuple[SnpHotspotRecord, ...]]:
    records_by_gene: dict[str, list[SnpHotspotRecord]] = {}

    for row in rows:
        gene = row.get("Hugo_Symbol", "")
        position = _parse_int(row.get("Amino_Acid_Position"))
        mutation_count = _parse_int(row.get("Mutation_Count"))
        alt_amino_acid, variant_count = _split_counted_value(row.get("Variant_Amino_Acid"))
        ref_amino_acid = row.get("ref") or _split_counted_value(row.get("Reference_Amino_Acid"))[0]

        if not gene or position is None or mutation_count is None:
            continue
        if not ref_amino_acid or not alt_amino_acid or variant_count is None:
            continue
        if len(ref_amino_acid) != 1 or len(alt_amino_acid) != 1:
            continue

        records_by_gene.setdefault(gene, []).append(
            SnpHotspotRecord(
                gene=gene,
                position=position,
                ref_amino_acid=ref_amino_acid,
                alt_amino_acid=alt_amino_acid,
                mutation_count=mutation_count,
                variant_count=variant_count,
            )
        )

    return {gene: tuple(records) for gene, records in records_by_gene.items()}


def _build_indel_index(rows: list[dict[str, str]]) -> dict[str, tuple[IndelHotspotRecord, ...]]:
    records_by_gene: dict[str, list[IndelHotspotRecord]] = {}

    for row in rows:
        gene = row.get("Hugo_Symbol", "")
        mutation_count = _parse_int(row.get("Mutation_Count"))
        variant_event, variant_count = _split_counted_value(row.get("Variant_Amino_Acid"))
        hotspot_position = row.get("Amino_Acid_Position", "")

        if not gene or mutation_count is None or not variant_event or variant_count is None:
            continue

        event_start, event_end = _parse_protein_event_bounds(variant_event)

        records_by_gene.setdefault(gene, []).append(
            IndelHotspotRecord(
                gene=gene,
                hotspot_position=hotspot_position,
                event=variant_event,
                event_start=event_start,
                event_end=event_end,
                mutation_count=mutation_count,
                variant_count=variant_count,
            )
        )

    return {gene: tuple(records) for gene, records in records_by_gene.items()}


def _split_counted_value(value: str | None) -> tuple[str | None, int | None]:
    if not value:
        return None, None

    parts = value.split(":", maxsplit=1)
    label = parts[0].strip()
    count = _parse_int(parts[1]) if len(parts) == 2 else None
    return (label or None), count


def _parse_int(value: str | None) -> int | None:
    if value is None:
        return None

    stripped = value.strip()
    if not stripped:
        return None

    try:
        return int(float(stripped))
    except ValueError:
        return None


def _match_transcript_consequence(
    gene_symbol: str,
    consequence: TranscriptConsequence,
    hotspot_index: HotspotIndex,
) -> EvidenceResult | None:
    if consequence.proteinEventType == "substitution":
        return _match_snv_consequence(gene_symbol, consequence, hotspot_index)

    if consequence.proteinEventType in INDIRECT_EVENT_TYPES:
        return _match_indel_consequence(gene_symbol, consequence, hotspot_index)

    return None


def _match_snv_consequence(
    gene_symbol: str,
    consequence: TranscriptConsequence,
    hotspot_index: HotspotIndex,
) -> EvidenceResult | None:
    if not consequence.proteinHgvs:
        return None

    match = re.match(r"^p\.([A-Z*])(\d+)([A-Z*=])$", consequence.proteinHgvs)
    if not match:
        return None

    ref_amino_acid, position_text, alt_amino_acid = match.groups()
    position = int(position_text)

    for record in hotspot_index.snv_by_gene.get(gene_symbol, ()):
        if (
            record.position == position
            and record.ref_amino_acid == ref_amino_acid
            and record.alt_amino_acid == alt_amino_acid
        ):
            return _build_hotspot_match_result(
                record.mutation_count,
                record.variant_count,
                consequence,
                {
                    "sheet": "SNV-hotspots",
                    "gene": gene_symbol,
                    "position": record.position,
                    "refAminoAcid": record.ref_amino_acid,
                    "altAminoAcid": record.alt_amino_acid,
                },
            )

    return None


def _match_indel_consequence(
    gene_symbol: str,
    consequence: TranscriptConsequence,
    hotspot_index: HotspotIndex,
) -> EvidenceResult | None:
    if not consequence.proteinHgvs:
        return None

    event_key = consequence.proteinHgvs.removeprefix("p.")
    if not event_key:
        return None

    for record in hotspot_index.indel_by_gene.get(gene_symbol, ()):
        if record.event != event_key:
            continue
        if not _positions_compatible(consequence, record):
            continue

        return _build_hotspot_match_result(
            record.mutation_count,
            record.variant_count,
            consequence,
            {
                "sheet": "INDEL-hotspots",
                "gene": gene_symbol,
                "position": record.hotspot_position,
                "event": record.event,
            },
        )

    return None


def _positions_compatible(
    consequence: TranscriptConsequence,
    record: IndelHotspotRecord,
) -> bool:
    consequence_start, consequence_end = _parse_consequence_event_bounds(consequence)

    if consequence_start is None or consequence_end is None:
        if consequence.proteinStart is None or consequence.proteinEnd is None:
            return True
        consequence_start = consequence.proteinStart
        consequence_end = consequence.proteinEnd

    if record.event_start is not None and record.event_end is not None:
        return (
            consequence_start == record.event_start
            and consequence_end == record.event_end
        )

    hotspot_start, hotspot_end = _parse_hotspot_position_bounds(record.hotspot_position)
    if hotspot_start is not None and hotspot_end is not None:
        return hotspot_start <= consequence_start and consequence_end <= hotspot_end

    return True


def _parse_consequence_event_bounds(
    consequence: TranscriptConsequence,
) -> tuple[int | None, int | None]:
    if not consequence.proteinHgvs:
        return None, None

    return _parse_protein_event_bounds(consequence.proteinHgvs.removeprefix("p."))


def _parse_protein_event_bounds(event: str) -> tuple[int | None, int | None]:
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


def _parse_hotspot_position_bounds(position: str) -> tuple[int | None, int | None]:
    if not position:
        return None, None

    range_match = re.match(r"^(\d+)-(\d+)$", position)
    if range_match:
        start_text, end_text = range_match.groups()
        return int(start_text), int(end_text)

    single_position = _parse_int(position)
    if single_position is not None:
        return single_position, single_position

    return None, None


def _build_hotspot_match_result(
    mutation_count: int,
    variant_count: int,
    consequence: TranscriptConsequence,
    matched_data: dict[str, object],
) -> EvidenceResult:
    enriched_matched_data = {
        **matched_data,
        "transcriptRefSeq": consequence.transcriptRefSeq,
        "proteinHgvs": consequence.proteinHgvs,
        "proteinEventType": consequence.proteinEventType,
        "proteinStart": consequence.proteinStart,
        "proteinEnd": consequence.proteinEnd,
        "mutationCount": mutation_count,
        "variantCount": variant_count,
    }

    if mutation_count >= 50 and variant_count >= 10:
        return EvidenceResult(
            score=4,
            evidenceCode="OS3",
            evidenceStatement=(
                "Located in Cancer Hotspots with at least 50 observed mutations "
                f"({mutation_count}) and at least 10 occurrences of the same protein event ({variant_count})."
            ),
            status="applied",
            source="cancerhotspots",
            matchedData=enriched_matched_data,
            dataAbsentReason=None,
        )

    if variant_count >= 10:
        return EvidenceResult(
            score=2,
            evidenceCode="OM3",
            evidenceStatement=(
                "Located in Cancer Hotspots with fewer than 50 total observed mutations "
                f"({mutation_count}), but at least 10 occurrences of the same protein event ({variant_count})."
            ),
            status="applied",
            source="cancerhotspots",
            matchedData=enriched_matched_data,
            dataAbsentReason=None,
        )

    return EvidenceResult(
        score=1,
        evidenceCode="OP3",
        evidenceStatement=(
            "Located in Cancer Hotspots, but the same protein event was observed fewer than 10 times "
            f"({variant_count})."
        ),
        status="applied",
        source="cancerhotspots",
        matchedData=enriched_matched_data,
        dataAbsentReason=None,
    )
