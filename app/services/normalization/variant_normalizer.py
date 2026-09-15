import re

import httpx

from app.models.normalized_variant import (
    AssemblyCoordinates,
    CoordinateRepresentations,
    GenomicHgvsRepresentations,
    NormalizationMetadata,
    NormalizedVariant,
    ProteinRepresentations,
    TranscriptHgvsRepresentations,
    VariantIdentifiers,
)


CLINGEN_ALLELE_REGISTRY_URL = "https://reg.clinicalgenome.org/allele"
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


class VariantNormalizationError(Exception):
    pass


def fetch_allele_registry_record(submitted_variant: str) -> dict:
    try:
        response = httpx.get(
            CLINGEN_ALLELE_REGISTRY_URL,
            params={"hgvs": submitted_variant},
            timeout=30.0,
        )
        response.raise_for_status()
    except httpx.HTTPError as exc:
        raise VariantNormalizationError(
            f"ClinGen normalization failed for '{submitted_variant}'."
        ) from exc

    payload = response.json()
    if not isinstance(payload, dict) or not payload:
        raise VariantNormalizationError(
            f"ClinGen returned no normalization result for '{submitted_variant}'."
        )

    return payload


def normalize_variant(submitted_variant: str) -> NormalizedVariant:
    record = fetch_allele_registry_record(submitted_variant)
    gene_symbol, gene_ncbi_id = _select_gene_metadata(record)

    genomic_hgvs, coordinates = _extract_genomic_representations(record)
    transcript_hgvs = _extract_transcript_representations(record)
    protein = _extract_protein_representations(record, gene_symbol)

    return NormalizedVariant(
        submitted_variant=submitted_variant,
        normalization=NormalizationMetadata(
            source="ClinGen Allele Registry",
            queried_variant=submitted_variant,
        ),
        identifiers=VariantIdentifiers(caid=_extract_caid(record)),
        geneSymbol=gene_symbol,
        geneNCBI_id=gene_ncbi_id,
        genomic_hgvs=genomic_hgvs,
        transcript_hgvs=transcript_hgvs,
        protein=protein,
        coordinates=coordinates,
    )


def _extract_caid(record: dict) -> str | None:
    allele_identifier = record.get("@id")
    if not isinstance(allele_identifier, str):
        return None
    return allele_identifier.rstrip("/").rsplit("/", maxsplit=1)[-1]


def _select_gene_metadata(record: dict) -> tuple[str | None, int | None]:
    transcript_alleles = record.get("transcriptAlleles") or []

    for transcript_allele in transcript_alleles:
        mane = transcript_allele.get("MANE") or {}
        if mane.get("maneStatus") == "MANE Select" and transcript_allele.get("geneSymbol"):
            return transcript_allele["geneSymbol"], transcript_allele.get("geneNCBI_id")

    for transcript_allele in transcript_alleles:
        gene_symbol = transcript_allele.get("geneSymbol")
        if gene_symbol:
            return gene_symbol, transcript_allele.get("geneNCBI_id")

    return None, None


def _extract_genomic_representations(
    record: dict,
) -> tuple[GenomicHgvsRepresentations, CoordinateRepresentations]:
    genomic_hgvs = GenomicHgvsRepresentations()
    coordinates = CoordinateRepresentations()

    for genomic_allele in record.get("genomicAlleles") or []:
        assembly = genomic_allele.get("referenceGenome")
        if assembly not in {"GRCh38", "GRCh37"}:
            continue

        hgvs_value = _select_preferred_hgvs(genomic_allele.get("hgvs") or [], ("NC_",))
        if assembly == "GRCh38":
            genomic_hgvs.GRCh38 = hgvs_value
            coordinates.GRCh38 = _extract_coordinates(genomic_allele, hgvs_value)
        if assembly == "GRCh37":
            genomic_hgvs.GRCh37 = hgvs_value
            coordinates.GRCh37 = _extract_coordinates(genomic_allele, hgvs_value)

    return genomic_hgvs, coordinates


def _extract_coordinates(genomic_allele: dict, hgvs_value: str | None) -> AssemblyCoordinates:
    coordinate_blocks = genomic_allele.get("coordinates") or []
    coordinate = coordinate_blocks[0] if coordinate_blocks else {}
    chromosome = genomic_allele.get("chromosome")
    refseq = hgvs_value.split(":", maxsplit=1)[0] if hgvs_value else None

    return AssemblyCoordinates(
        refseq=refseq,
        chrom=_normalize_chrom_label(chromosome),
        chrom_num=_normalize_chrom_num(chromosome),
        pos=coordinate.get("end"),
        ref=coordinate.get("referenceAllele"),
        alt=coordinate.get("allele"),
    )


def _normalize_chrom_label(chromosome: str | None) -> str | None:
    if chromosome == "MT":
        return "chrM"
    if chromosome:
        return f"chr{chromosome}"

    return None


def _normalize_chrom_num(chromosome: str | None) -> str | None:
    if chromosome == "X":
        return "23"
    if chromosome == "Y":
        return "24"
    if chromosome == "MT":
        return "M"

    return chromosome


def _extract_transcript_representations(record: dict) -> TranscriptHgvsRepresentations:
    transcript_hgvs = TranscriptHgvsRepresentations()
    transcript_alleles = record.get("transcriptAlleles") or []
    first_nm_hgvs = _select_first_nm_hgvs(transcript_alleles)

    for transcript_allele in transcript_alleles:
        mane = transcript_allele.get("MANE") or {}
        if mane.get("maneStatus") == "MANE Select":
            refseq_hgvs = (((mane.get("nucleotide") or {}).get("RefSeq") or {}).get("hgvs"))
            if refseq_hgvs and refseq_hgvs.startswith("NM_"):
                transcript_hgvs.mane_select_b38 = refseq_hgvs
                transcript_hgvs.mane_select_b38_source = "clingen"
                break

    for transcript_allele in transcript_alleles:
        genome_alignments = transcript_allele.get("genomeAlignments") or []
        has_grch37_alignment = any(
            alignment.get("referenceGenome") == "GRCh37"
            for alignment in genome_alignments
        )
        if not has_grch37_alignment:
            continue

        for hgvs_value in transcript_allele.get("hgvs") or []:
            if hgvs_value.startswith("NM_"):
                transcript_hgvs.canonical_b37 = hgvs_value
                break

        if transcript_hgvs.canonical_b37:
            break

    transcript_hgvs.representative_transcript_hgvs = (
        transcript_hgvs.mane_select_b38
        or transcript_hgvs.canonical_b37
        or first_nm_hgvs
    )

    return transcript_hgvs


def _select_first_nm_hgvs(transcript_alleles: list[dict]) -> str | None:
    for transcript_allele in transcript_alleles:
        for hgvs_value in transcript_allele.get("hgvs") or []:
            if hgvs_value.startswith("NM_"):
                return hgvs_value

    return None


def _extract_protein_representations(record: dict, gene: str | None) -> ProteinRepresentations:
    protein = ProteinRepresentations()
    transcript_alleles = record.get("transcriptAlleles") or []

    protein_hgvs_full = None
    for transcript_allele in transcript_alleles:
        mane = transcript_allele.get("MANE") or {}
        if mane.get("maneStatus") == "MANE Select":
            protein_hgvs_full = (((mane.get("protein") or {}).get("RefSeq") or {}).get("hgvs"))
            if protein_hgvs_full and protein_hgvs_full.startswith("NP_"):
                break
            protein_hgvs_full = None

    if not protein_hgvs_full:
        for transcript_allele in transcript_alleles:
            protein_effect = transcript_allele.get("proteinEffect") or {}
            candidate_hgvs = protein_effect.get("hgvsWellDefined") or protein_effect.get("hgvs")
            if candidate_hgvs and candidate_hgvs.startswith("NP_"):
                protein_hgvs_full = candidate_hgvs
                break

    if protein_hgvs_full:
        protein.hgvs_protein_full = protein_hgvs_full
        protein.hgvs_3letter = _strip_accession_prefix(protein_hgvs_full)
        protein.hgvs_1letter = _convert_three_letter_protein_hgvs(protein.hgvs_3letter)
        protein.short_name = (
            protein.hgvs_1letter.removeprefix("p.") if protein.hgvs_1letter else None
        )
        protein.np_accession = (
            protein_hgvs_full.split(":", maxsplit=1)[0]
            if protein_hgvs_full.startswith("NP_")
            else None
        )

    if gene and protein.short_name:
        protein.civic_profile_name = f"{gene} {protein.short_name}"

    return protein


def _strip_accession_prefix(hgvs_value: str) -> str:
    return hgvs_value.split(":", maxsplit=1)[-1]


def _convert_three_letter_protein_hgvs(hgvs_value: str | None) -> str | None:
    if not hgvs_value:
        return None

    def replace_match(match: re.Match[str]) -> str:
        return AMINO_ACID_THREE_TO_ONE.get(match.group(0), match.group(0))

    return re.sub(
        r"Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter",
        replace_match,
        hgvs_value,
    )


def _select_preferred_hgvs(hgvs_values: list[str], prefixes: tuple[str, ...]) -> str | None:
    for hgvs_value in hgvs_values:
        if hgvs_value.startswith(prefixes):
            return hgvs_value
    return None
