"""Harvest ClinGen CSPEC PM1/OM1 text into an OM1 working CSV.

This utility is intentionally conservative.

ClinGen CSPEC records are the source of truth for v1 OM1 curation, but the
released JSON does not consistently expose functional domain boundaries as clean
structured fields. In practice, many specifications encode residue ranges in the
PM1 evidence-strength description text. This script downloads released CSPEC
records for a target gene set, enriches each row with the current MANE Select
RefSeq transcript and protein accessions, extracts any obvious residue intervals
or codon lists, and writes a CSV in the planned OM1 row shape.

The generated CSV is not intended to be consumed directly by scoring code. It is
an intermediate curation sheet that follows the planned OM1 row shape while
embedding CSPEC provenance and raw rule text into the notes column. Every target
gene is represented in the output. The `rowStatus` column distinguishes rows
that are operational today from rows that should currently be ignored, such as
genes with no released ClinGen CSPEC specification, genes with only non-released
ClinGen specifications, or genes whose published CSPEC text lacks
machine-readable domain boundaries.

Current `rowStatus` values:

- `ready`: ClinGen text exposed a usable interval or codon list.
- `missing_boundaries`: a released ClinGen rule exists, but no machine-readable
    interval could be extracted from the published CSPEC text.
- `clingen_spec_not_released`: ClinGen has a registry record for the gene, but
    the detailed specification is not in an allowed release state for this build.
- `no_clingen_spec_found`: no ClinGen CSPEC record was found for the gene in the
    registry.

Usage examples:

    python3 -m app.services.utilities.clingen_om1_seed_builder
    python3 -m app.services.utilities.clingen_om1_seed_builder --gene BRAF --gene KRAS
    python3 -m app.services.utilities.clingen_om1_seed_builder --out data/om1_clingen_domains_seed.csv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import ssl
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


CSPEC_INDEX_URL = "https://cspec.genome.network/cspec/api/svis"
MANE_SUMMARY_DIR_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
DEFAULT_TIMEOUT_SECONDS = 30

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_GENES_FILE = REPO_ROOT / "data" / "clinmave" / "genes.txt"
DEFAULT_OUTPUT_FILE = REPO_ROOT / "data" / "om1_clingen_domains_seed.csv"
RELEASED_STATUSES = {"Released"}
RULE_LABELS = {"PM1", "OM1"}

AA_RANGE_PATTERN = re.compile(
    r"(?P<name>[^,;()\[]+?)\s*\[AA\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)\]",
    flags=re.IGNORECASE,
)
ACCESSION_PATTERN = re.compile(r"\b([A-Z]{2}_[\d .]+(?:\.\d+)?)\b")
RANGE_LIST_PATTERN = re.compile(r"(?P<start>\d+)\s*-\s*(?P<end>\d+)")
CODON_LIST_PATTERN = re.compile(
    r"codons?(?:[^:]{0,80})?:\s*(?P<values>\d+(?:\s*,\s*\d+)*)",
    flags=re.IGNORECASE,
)
MOTIF_LIST_PATTERN = re.compile(
    r"residues?\s+in\s+(?P<label>[^:]+):\s*(?P<values>\d+\s*-\s*\d+(?:\s*,\s*\d+\s*-\s*\d+)*)",
    flags=re.IGNORECASE,
)
EXCLUSION_PATTERN = re.compile(
    r"exclude(?:d|s|ing)?\s+(?:residues?|codons?)\s*[:]?\s*(?P<values>[^.;]+)",
    flags=re.IGNORECASE,
)


@dataclass(slots=True)
class SeedRow:
    gene_symbol: str
    criteria_label: str
    evidence_strength: str
    spec_id: str
    spec_api_url: str
    spec_ui_url: str
    spec_status: str
    spec_version: str
    source_transcript_accession: str
    source_protein_accession: str
    domain_name: str
    start_residue: str
    end_residue: str
    excluded_residues: str
    row_status: str
    extraction_status: str
    raw_description: str
    notes: str


@dataclass(slots=True)
class ManeReference:
    gene_symbol: str
    refseq_nuc: str
    refseq_prot: str
    mane_status: str
    summary_version: str


class CSpecClient:
    def __init__(self, timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS) -> None:
        self.timeout_seconds = timeout_seconds
        self.ssl_context = ssl.create_default_context()
        self.headers = {"User-Agent": "Mozilla/5.0"}

    def get_json(self, url: str) -> Any:
        request = urllib.request.Request(url, headers=self.headers)
        with urllib.request.urlopen(
            request,
            timeout=self.timeout_seconds,
            context=self.ssl_context,
        ) as response:
            return json.load(response)

    def get_text(self, url: str) -> str:
        request = urllib.request.Request(url, headers=self.headers)
        with urllib.request.urlopen(
            request,
            timeout=self.timeout_seconds,
            context=self.ssl_context,
        ) as response:
            return response.read().decode("utf-8", "replace")

    def get_bytes(self, url: str) -> bytes:
        request = urllib.request.Request(url, headers=self.headers)
        with urllib.request.urlopen(
            request,
            timeout=self.timeout_seconds,
            context=self.ssl_context,
        ) as response:
            return response.read()


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a ClinGen-only OM1 working spreadsheet from CSPEC and MANE Select.",
    )
    parser.add_argument(
        "--gene",
        action="append",
        dest="genes",
        default=[],
        help="Gene symbol to include. Repeat to include multiple genes.",
    )
    parser.add_argument(
        "--genes-file",
        type=Path,
        default=DEFAULT_GENES_FILE,
        help=f"Optional file with one gene symbol per line. Defaults to {DEFAULT_GENES_FILE}.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT_FILE,
        help=f"Output CSV path in the planned OM1 row shape. Defaults to {DEFAULT_OUTPUT_FILE}.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=DEFAULT_TIMEOUT_SECONDS,
        help=f"HTTP timeout in seconds. Defaults to {DEFAULT_TIMEOUT_SECONDS}.",
    )
    parser.add_argument(
        "--include-status",
        action="append",
        dest="statuses",
        default=[],
        help="Additional CSPEC statuses to include beyond Released.",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def find_current_mane_summary_filename(client: CSpecClient) -> str:
    listing = client.get_text(MANE_SUMMARY_DIR_URL)
    match = re.search(r'href="(MANE\.GRCh38\.v[^"/]+\.summary\.txt\.gz)"', listing)
    if not match:
        raise ValueError("Could not locate current MANE summary file")
    return match.group(1)


def parse_mane_summary_version(filename: str) -> str:
    match = re.search(r"MANE\.GRCh38\.(v[0-9.]+)\.summary\.txt\.gz", filename)
    return match.group(1) if match else ""


def load_mane_reference_map(client: CSpecClient) -> dict[str, ManeReference]:
    filename = find_current_mane_summary_filename(client)
    summary_version = parse_mane_summary_version(filename)
    payload = gzip.decompress(client.get_bytes(f"{MANE_SUMMARY_DIR_URL}{filename}")).decode(
        "utf-8",
        "replace",
    )
    reader = csv.DictReader(io.StringIO(payload), delimiter="\t")

    by_gene: dict[str, list[dict[str, str]]] = {}
    for row in reader:
        symbol = (row.get("symbol") or "").strip()
        if not symbol:
            continue
        by_gene.setdefault(symbol, []).append(row)

    mane_by_gene: dict[str, ManeReference] = {}
    for symbol, rows in by_gene.items():
        selected_row = next((row for row in rows if row.get("MANE_status") == "MANE Select"), None)
        if selected_row is None:
            continue
        mane_by_gene[symbol] = ManeReference(
            gene_symbol=symbol,
            refseq_nuc=(selected_row.get("RefSeq_nuc") or "").strip(),
            refseq_prot=(selected_row.get("RefSeq_prot") or "").strip(),
            mane_status=(selected_row.get("MANE_status") or "").strip(),
            summary_version=summary_version,
        )
    return mane_by_gene


def read_gene_file(path: Path | None) -> list[str]:
    if path is None or not path.exists():
        return []
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def load_target_genes(args: argparse.Namespace) -> list[str]:
    genes = [gene.strip() for gene in args.genes if gene and gene.strip()]
    genes.extend(read_gene_file(args.genes_file))
    return sorted(set(genes))


def unescape_cspec_text(text: str) -> str:
    return (
        text.replace(r"\_", "_")
        .replace(r"\[", "[")
        .replace(r"\]", "]")
        .replace(r"\<", "<")
        .replace(r"\>", ">")
    )


def normalize_accession(value: str) -> str:
    return value.replace(" ", "")


def extract_accessions(text: str) -> tuple[str, str]:
    transcript = ""
    protein = ""
    for match in ACCESSION_PATTERN.findall(text):
        accession = normalize_accession(match)
        if accession.startswith("NM_") and not transcript:
            transcript = accession
        if accession.startswith("NP_") and not protein:
            protein = accession
    return transcript, protein


def extract_excluded_residues(text: str) -> str:
    match = EXCLUSION_PATTERN.search(text)
    if not match:
        return ""
    return " ".join(match.group("values").split())


def parse_named_ranges(text: str) -> list[tuple[str, int, int]]:
    results: list[tuple[str, int, int]] = []
    for match in AA_RANGE_PATTERN.finditer(text):
        domain_name = " ".join(match.group("name").split()).strip(" ,;:")
        start = int(match.group("start"))
        end = int(match.group("end"))
        results.append((domain_name, start, end))
    return results


def parse_motif_ranges(text: str) -> list[tuple[str, int, int]]:
    match = MOTIF_LIST_PATTERN.search(text)
    if not match:
        return []
    label = " ".join(match.group("label").split()).strip(" ,;:")
    results: list[tuple[str, int, int]] = []
    for range_match in RANGE_LIST_PATTERN.finditer(match.group("values")):
        start = int(range_match.group("start"))
        end = int(range_match.group("end"))
        results.append((label, start, end))
    return results


def parse_codon_list(text: str) -> list[int]:
    match = CODON_LIST_PATTERN.search(text)
    if not match:
        return []
    return [int(value.strip()) for value in match.group("values").split(",") if value.strip()]


def build_seed_rows(
    gene_symbol: str,
    criteria_label: str,
    evidence_strength: str,
    description: str,
    spec_id: str,
    spec_api_url: str,
    spec_ui_url: str,
    spec_status: str,
    spec_version: str,
) -> list[SeedRow]:
    normalized_text = " ".join(unescape_cspec_text(description).split())
    transcript, protein = extract_accessions(normalized_text)
    excluded_residues = extract_excluded_residues(normalized_text)

    rows: list[SeedRow] = []
    for domain_name, start, end in parse_named_ranges(normalized_text):
        rows.append(
            SeedRow(
                gene_symbol=gene_symbol,
                criteria_label=criteria_label,
                evidence_strength=evidence_strength,
                spec_id=spec_id,
                spec_api_url=spec_api_url,
                spec_ui_url=spec_ui_url,
                spec_status=spec_status,
                spec_version=spec_version,
                source_transcript_accession=transcript,
                source_protein_accession=protein,
                domain_name=domain_name,
                start_residue=str(start),
                end_residue=str(end),
                excluded_residues=excluded_residues,
                row_status="ready",
                extraction_status="parsed_range",
                raw_description=normalized_text,
                notes="Parsed named AA range from CSPEC text.",
            )
        )
    if rows:
        return rows

    motif_rows = []
    for domain_name, start, end in parse_motif_ranges(normalized_text):
        motif_rows.append(
            SeedRow(
                gene_symbol=gene_symbol,
                criteria_label=criteria_label,
                evidence_strength=evidence_strength,
                spec_id=spec_id,
                spec_api_url=spec_api_url,
                spec_ui_url=spec_ui_url,
                spec_status=spec_status,
                spec_version=spec_version,
                source_transcript_accession=transcript,
                source_protein_accession=protein,
                domain_name=domain_name,
                start_residue=str(start),
                end_residue=str(end),
                excluded_residues=excluded_residues,
                row_status="ready",
                extraction_status="parsed_motif_range",
                raw_description=normalized_text,
                notes="Parsed unlabeled residue range list from CSPEC text.",
            )
        )
    if motif_rows:
        return motif_rows

    codon_rows = []
    for residue in parse_codon_list(normalized_text):
        codon_rows.append(
            SeedRow(
                gene_symbol=gene_symbol,
                criteria_label=criteria_label,
                evidence_strength=evidence_strength,
                spec_id=spec_id,
                spec_api_url=spec_api_url,
                spec_ui_url=spec_ui_url,
                spec_status=spec_status,
                spec_version=spec_version,
                source_transcript_accession=transcript,
                source_protein_accession=protein,
                domain_name="Curated codon hotspot",
                start_residue=str(residue),
                end_residue=str(residue),
                excluded_residues=excluded_residues,
                row_status="ready",
                extraction_status="parsed_codon",
                raw_description=normalized_text,
                notes="Parsed codon list from CSPEC text.",
            )
        )
    if codon_rows:
        return codon_rows

    return [
        SeedRow(
            gene_symbol=gene_symbol,
            criteria_label=criteria_label,
            evidence_strength=evidence_strength,
            spec_id=spec_id,
            spec_api_url=spec_api_url,
            spec_ui_url=spec_ui_url,
            spec_status=spec_status,
            spec_version=spec_version,
            source_transcript_accession=transcript,
            source_protein_accession=protein,
            domain_name="",
            start_residue="",
            end_residue="",
            excluded_residues=excluded_residues,
            row_status="missing_boundaries",
            extraction_status="manual_review_required",
            raw_description=normalized_text,
            notes="No explicit residue intervals could be parsed automatically.",
        )
    ]


def build_inventory_placeholder_row(
    gene_symbol: str,
    row_status: str,
    notes: str,
    *,
    spec_id: str = "",
    spec_api_url: str = "",
    spec_ui_url: str = "",
    spec_status: str = "",
    spec_version: str = "",
) -> SeedRow:
    return SeedRow(
        gene_symbol=gene_symbol,
        criteria_label="",
        evidence_strength="",
        spec_id=spec_id,
        spec_api_url=spec_api_url,
        spec_ui_url=spec_ui_url,
        spec_status=spec_status,
        spec_version=spec_version,
        source_transcript_accession="",
        source_protein_accession="",
        domain_name="",
        start_residue="",
        end_residue="",
        excluded_residues="",
        row_status=row_status,
        extraction_status="",
        raw_description="",
        notes=notes,
    )


def collect_specifications(
    client: CSpecClient,
    target_genes: set[str],
    allowed_statuses: set[str],
) -> list[dict[str, Any]]:
    payload = client.get_json(CSPEC_INDEX_URL)
    if not isinstance(payload, dict) or not isinstance(payload.get("data"), list):
        raise ValueError("Unexpected CSPEC index response shape")

    matches: list[dict[str, Any]] = []
    for entry in payload["data"]:
        status = entry.get("status") or entry.get("currentStatus") or ""
        if status not in allowed_statuses:
            continue
        genes = []
        for ruleset in entry.get("ruleSets", []):
            for gene in ruleset.get("genes", []):
                label = gene.get("label")
                if label:
                    genes.append(label)
        if target_genes and not target_genes.intersection(genes):
            continue
        matches.append(entry)
    return matches


def collect_all_target_gene_spec_statuses(
    client: CSpecClient,
    target_genes: set[str],
) -> dict[str, list[tuple[str, str]]]:
    payload = client.get_json(CSPEC_INDEX_URL)
    if not isinstance(payload, dict) or not isinstance(payload.get("data"), list):
        raise ValueError("Unexpected CSPEC index response shape")

    statuses_by_gene: dict[str, list[tuple[str, str]]] = {}
    for entry in payload["data"]:
        spec_api_url = entry.get("@id")
        if not spec_api_url:
            continue
        detailed_spec = client.get_json(spec_api_url)
        detailed_status = (
            detailed_spec.get("currentStatus")
            or detailed_spec.get("cspecStatus")
            or detailed_spec.get("status")
            or ""
        )
        spec_id = spec_api_url.rsplit("/", 1)[-1]
        for gene_symbol in collect_spec_genes(detailed_spec, target_genes):
            statuses_by_gene.setdefault(gene_symbol, []).append((spec_id, detailed_status))
    return statuses_by_gene


def iter_rule_rows(spec: dict[str, Any], target_genes: set[str]) -> list[SeedRow]:
    spec_api_url = spec["@id"]
    spec_id = spec_api_url.rsplit("/", 1)[-1]
    spec_ui_url = f"https://cspec.genome.network/cspec/ui/svi/doc/{spec_id}"
    spec_status = spec.get("currentStatus") or spec.get("cspecStatus") or spec.get("status") or ""
    spec_version = str(spec.get("version") or "")

    rows: list[SeedRow] = []
    for ruleset in spec.get("ruleSets", []):
        ruleset_genes = {
            gene.get("label")
            for gene in ruleset.get("genes", [])
            if gene.get("label")
        }
        matched_genes = sorted(ruleset_genes & target_genes) if target_genes else sorted(ruleset_genes)
        if not matched_genes:
            continue

        for criteria_code in ruleset.get("criteriaCodes", []):
            criteria_label = criteria_code.get("label")
            if criteria_label not in RULE_LABELS:
                continue

            applicable_strengths = [
                strength
                for strength in criteria_code.get("evidenceStrengths", [])
                if str(strength.get("applicability", "")).lower() == "applicable"
                and strength.get("description")
            ]

            descriptions: list[tuple[str, str]] = []
            if applicable_strengths:
                descriptions.extend(
                    (str(strength.get("label") or ""), str(strength.get("description") or ""))
                    for strength in applicable_strengths
                )
            elif criteria_code.get("description"):
                descriptions.append(("criterion", str(criteria_code.get("description") or "")))

            for gene_symbol in matched_genes:
                for evidence_strength, description in descriptions:
                    rows.extend(
                        build_seed_rows(
                            gene_symbol=gene_symbol,
                            criteria_label=criteria_label,
                            evidence_strength=evidence_strength,
                            description=description,
                            spec_id=spec_id,
                            spec_api_url=spec_api_url,
                            spec_ui_url=spec_ui_url,
                            spec_status=spec_status,
                            spec_version=spec_version,
                        )
                    )
    return rows


def collect_spec_genes(spec: dict[str, Any], target_genes: set[str]) -> set[str]:
    matched_genes: set[str] = set()
    for ruleset in spec.get("ruleSets", []):
        for gene in ruleset.get("genes", []):
            label = gene.get("label")
            if label and (not target_genes or label in target_genes):
                matched_genes.add(label)
    return matched_genes


def write_seed_csv(path: Path, rows: list[SeedRow], mane_by_gene: dict[str, ManeReference]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            mane_reference = mane_by_gene.get(row.gene_symbol)
            mane_transcript = mane_reference.refseq_nuc if mane_reference else ""
            mane_protein = mane_reference.refseq_prot if mane_reference else ""
            explicit_reference_notes = []
            if row.source_transcript_accession:
                explicit_reference_notes.append(f"clingenTranscript={row.source_transcript_accession}")
            if row.source_protein_accession:
                explicit_reference_notes.append(f"clingenProtein={row.source_protein_accession}")
            notes = " | ".join(
                part
                for part in [
                    row.notes,
                    f"rowStatus={row.row_status}",
                    f"criteria={row.criteria_label}",
                    f"strength={row.evidence_strength}",
                    f"status={row.spec_status}",
                    f"specId={row.spec_id}",
                    f"extraction={row.extraction_status}",
                    f"maneStatus={mane_reference.mane_status}" if mane_reference else "maneStatus=missing",
                    f"maneVersion={mane_reference.summary_version}" if mane_reference else "",
                    f"specUrl={row.spec_ui_url}",
                    *explicit_reference_notes,
                    f"rawDescription={row.raw_description}",
                ]
                if part
            )
            writer.writerow(
                {
                    "geneSymbol": row.gene_symbol,
                    "maneTranscript": mane_transcript,
                    "maneProtein": mane_protein,
                    "domainName": row.domain_name,
                    "startResidue": row.start_residue,
                    "endResidue": row.end_residue,
                    "excludedResidues": row.excluded_residues,
                    "rowStatus": row.row_status,
                    "source": "ClinGen CSPEC",
                    "sourceVersion": row.spec_version,
                    "notes": notes,
                }
            )


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    target_genes = set(load_target_genes(args))
    allowed_statuses = set(RELEASED_STATUSES)
    allowed_statuses.update(status.strip() for status in args.statuses if status.strip())

    client = CSpecClient(timeout_seconds=args.timeout_seconds)
    matches = collect_specifications(client, target_genes, allowed_statuses)
    all_spec_statuses_by_gene = collect_all_target_gene_spec_statuses(client, target_genes)
    mane_by_gene = load_mane_reference_map(client)

    rows: list[SeedRow] = []
    genes_with_released_specs: set[str] = set()
    genes_with_rows: set[str] = set()
    for match in matches:
        detailed_spec = client.get_json(match["@id"])
        detailed_status = (
            detailed_spec.get("currentStatus")
            or detailed_spec.get("cspecStatus")
            or detailed_spec.get("status")
            or ""
        )
        if detailed_status not in allowed_statuses:
            continue

        matched_genes = collect_spec_genes(detailed_spec, target_genes)
        if not matched_genes:
            continue
        genes_with_released_specs.update(matched_genes)

        spec_rows = iter_rule_rows(detailed_spec, target_genes)
        if spec_rows:
            rows.extend(spec_rows)
            genes_with_rows.update(row.gene_symbol for row in spec_rows)
            continue

        spec_api_url = detailed_spec["@id"]
        spec_id = spec_api_url.rsplit("/", 1)[-1]
        spec_ui_url = f"https://cspec.genome.network/cspec/ui/svi/doc/{spec_id}"
        spec_version = str(detailed_spec.get("version") or "")
        for gene_symbol in sorted(matched_genes):
            rows.append(
                build_inventory_placeholder_row(
                    gene_symbol,
                    "no_applicable_rule",
                    "Released ClinGen CSPEC found, but no applicable PM1 or "
                    "OM1 rule text was available for extraction.",
                    spec_id=spec_id,
                    spec_api_url=spec_api_url,
                    spec_ui_url=spec_ui_url,
                    spec_status=detailed_status,
                    spec_version=spec_version,
                )
            )
            genes_with_rows.add(gene_symbol)

    for gene_symbol in sorted(target_genes - genes_with_released_specs):
        non_released_specs = all_spec_statuses_by_gene.get(gene_symbol, [])
        if non_released_specs:
            spec_summaries = ", ".join(
                f"{spec_id}:{status or 'unknown'}"
                for spec_id, status in sorted(non_released_specs)
            )
            rows.append(
                build_inventory_placeholder_row(
                    gene_symbol,
                    "clingen_spec_not_released",
                    "ClinGen CSPEC record exists, but no specification is "
                    "currently in an allowed release state for this build "
                    f"policy. Non-released specs: {spec_summaries}.",
                )
            )
            genes_with_rows.add(gene_symbol)
            continue
        rows.append(
            build_inventory_placeholder_row(
                gene_symbol,
                "no_clingen_spec_found",
                "No ClinGen CSPEC specification was found for this gene in the registry.",
            )
        )
        genes_with_rows.add(gene_symbol)

    for gene_symbol in sorted(genes_with_released_specs - genes_with_rows):
        rows.append(
            build_inventory_placeholder_row(
                gene_symbol,
                "no_inventory_row",
                "Released ClinGen CSPEC was found, but no inventory row was emitted. This indicates an extraction gap.",
            )
        )

    rows.sort(
        key=lambda row: (
            row.gene_symbol,
            row.row_status,
            row.criteria_label,
            row.extraction_status,
            int(row.start_residue) if row.start_residue.isdigit() else 999999,
            row.domain_name,
        )
    )
    write_seed_csv(args.out, rows, mane_by_gene)

    print(f"wrote {len(rows)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
