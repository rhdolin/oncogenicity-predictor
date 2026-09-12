from typing import Literal

from pydantic import BaseModel, Field


class NormalizationMetadata(BaseModel):
    source: str
    queried_variant: str


class VariantIdentifiers(BaseModel):
    caid: str | None = None


class GenomicHgvsRepresentations(BaseModel):
    GRCh38: str | None = None
    GRCh37: str | None = None


class TranscriptHgvsRepresentations(BaseModel):
    mane_select_b38: str | None = None
    mane_select_b38_source: Literal["clingen", "vep"] | None = None
    canonical_b37: str | None = None
    representative_transcript_hgvs: str | None = None


class ProteinRepresentations(BaseModel):
    np_accession: str | None = None
    hgvs_protein_full: str | None = None
    hgvs_3letter: str | None = None
    hgvs_1letter: str | None = None
    short_name: str | None = None
    civic_profile_name: str | None = None


class AssemblyCoordinates(BaseModel):
    refseq: str | None = None
    chrom: str | None = None
    chrom_num: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None


class CoordinateRepresentations(BaseModel):
    GRCh38: AssemblyCoordinates | None = None
    GRCh37: AssemblyCoordinates | None = None


class NormalizedVariant(BaseModel):
    submitted_variant: str
    normalization: NormalizationMetadata
    identifiers: VariantIdentifiers = Field(default_factory=VariantIdentifiers)
    geneSymbol: str | None = None
    geneNCBI_id: int | None = None
    genomic_hgvs: GenomicHgvsRepresentations = Field(default_factory=GenomicHgvsRepresentations)
    transcript_hgvs: TranscriptHgvsRepresentations = Field(default_factory=TranscriptHgvsRepresentations)
    protein: ProteinRepresentations = Field(default_factory=ProteinRepresentations)
    coordinates: CoordinateRepresentations = Field(default_factory=CoordinateRepresentations)


class NormalizedVariantBatchResponse(BaseModel):
    normalized_variants: list[NormalizedVariant]
