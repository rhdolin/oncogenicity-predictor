from typing import Literal

from pydantic import BaseModel, Field

from app.models.normalized_variant import NormalizedVariant


class PopulationSummary(BaseModel):
    maxSubpopulationAf: float | None = None
    maxSubpopulationLabel: str | None = None
    maxOverallAf: float | None = None
    maxOverallLabel: str | None = None


class TranscriptConsequence(BaseModel):
    transcriptRefSeq: str
    consequenceTerms: list[str] = Field(default_factory=list)
    proteinStart: int | None = None
    proteinEnd: int | None = None
    aminoAcids: str | None = None
    isManeSelect: bool | None = None


class BasicAnnotation(BaseModel):
    mostSevereConsequence: str | None = None
    transcriptConsequences: list[TranscriptConsequence] = Field(default_factory=list)
    population: PopulationSummary = Field(default_factory=PopulationSummary)


class CaddAnnotation(BaseModel):
    phred: float | None = None
    raw: float | None = None


class ComputationalAnnotation(BaseModel):
    cadd: CaddAnnotation = Field(default_factory=CaddAnnotation)
    phyloP100wayVertebrate: float | None = None


class AnnotationError(BaseModel):
    source: str
    message: str
    attemptedQueries: list[str] = Field(default_factory=list)


class AnnotatedVariant(BaseModel):
    normalizedVariant: NormalizedVariant
    annotationStatus: Literal["complete", "failed"]
    annotationError: AnnotationError | None = None
    basicAnnotation: BasicAnnotation | None = None
    computationalAnnotation: ComputationalAnnotation | None = None


class AnnotatedVariantBatchResponse(BaseModel):
    annotated_variants: list[AnnotatedVariant]
