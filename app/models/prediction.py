from typing import Any, Literal

from pydantic import BaseModel, Field


EvidenceStatus = Literal["applied", "not_available"]


class EvidenceResult(BaseModel):
    score: int
    evidenceCode: str | None = None
    evidenceStatement: str
    status: EvidenceStatus
    source: str | None = None
    matchedData: dict[str, Any] | None = None
    dataAbsentReason: str | None = None


class OncogenicityEvidence(BaseModel):
    population: EvidenceResult
    computational: EvidenceResult
    hotspots: EvidenceResult
    predictive: EvidenceResult
    functional: EvidenceResult


class OncogenicityPredictionSummary(BaseModel):
    overallScore: int
    overallClassification: str | None = None
    oncogenicityEvidence: OncogenicityEvidence


class Coding(BaseModel):
    system: str | None = None
    code: str | None = None
    display: str | None = None


class CodeableConcept(BaseModel):
    coding: list[Coding] = Field(default_factory=list)
    text: str | None = None


class Extension(BaseModel):
    url: str
    valueString: str | None = None


class ObservationComponent(BaseModel):
    code: CodeableConcept
    valueInteger: int | None = None
    interpretation: list[CodeableConcept] = Field(default_factory=list)
    dataAbsentReason: CodeableConcept | None = None


class OncogenicityObservation(BaseModel):
    resourceType: Literal["Observation"] = "Observation"
    status: Literal["final"] = "final"
    issued: str
    code: CodeableConcept
    valueInteger: int
    interpretation: list[CodeableConcept] = Field(default_factory=list)
    extension: list[Extension] = Field(default_factory=list)
    component: list[ObservationComponent] = Field(default_factory=list)


class OncogenicityPredictionBatchResponse(BaseModel):
    observations: list[OncogenicityObservation]
