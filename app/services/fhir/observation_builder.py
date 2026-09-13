from datetime import datetime, timezone

from app.models.prediction import (
    CodeableConcept,
    Coding,
    EvidenceResult,
    ObservationComponent,
    OncogenicityObservation,
    OncogenicityPredictionSummary,
    Reference,
)


TEMP_CODE_SYSTEM = "https://oncogenicity-predictor.example/fhir/CodeSystem/temp-codes"
DATA_ABSENT_REASON_SYSTEM = "http://terminology.hl7.org/CodeSystem/data-absent-reason"


def _build_concept(
    code: str | None = None,
    display: str | None = None,
    text: str | None = None,
) -> CodeableConcept:
    coding = []
    if code is not None or display is not None:
        coding.append(Coding(system=TEMP_CODE_SYSTEM, code=code, display=display))
    return CodeableConcept(coding=coding, text=text)


def _build_pipeline_component(
    pipeline_name: str,
    evidence: EvidenceResult,
) -> ObservationComponent:
    component = ObservationComponent(
        code=_build_concept(
            code=f"{pipeline_name}-evidence",
            display=f"{pipeline_name.title()} evidence",
            text=f"{pipeline_name.title()} evidence",
        )
    )

    if evidence.status == "not_available":
        component.dataAbsentReason = CodeableConcept(
            coding=[
                Coding(
                    system=DATA_ABSENT_REASON_SYSTEM,
                    code=evidence.dataAbsentReason,
                    display=evidence.dataAbsentReason,
                )
            ],
            text=evidence.evidenceStatement,
        )
        return component

    component.valueInteger = evidence.score
    component.interpretation = [
        _build_concept(
            code=evidence.evidenceCode,
            display=evidence.evidenceCode,
            text=evidence.evidenceStatement,
        )
    ]
    return component


def build_oncogenicity_observation(
    summary: OncogenicityPredictionSummary,
    submitted_variant: str,
) -> OncogenicityObservation:
    interpretation: list[CodeableConcept] = []
    if summary.overallClassification is not None:
        interpretation = [
            _build_concept(
                code=summary.overallClassification,
                display=summary.overallClassification,
                text=summary.overallClassification,
            )
        ]

    return OncogenicityObservation(
        issued=datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
            "+00:00",
            "Z",
        ),
        code=_build_concept(
            code="oncogenicity-prediction",
            display="Oncogenicity prediction",
            text="Oncogenicity prediction",
        ),
        valueInteger=summary.overallScore,
        interpretation=interpretation,
        derivedFrom=[Reference(reference=submitted_variant, display=submitted_variant)],
        component=[
            _build_pipeline_component("population", summary.oncogenicityEvidence.population),
            _build_pipeline_component(
                "computational",
                summary.oncogenicityEvidence.computational,
            ),
            _build_pipeline_component("hotspots", summary.oncogenicityEvidence.hotspots),
            _build_pipeline_component(
                "predictive",
                summary.oncogenicityEvidence.predictive,
            ),
            _build_pipeline_component("functional", summary.oncogenicityEvidence.functional),
        ],
    )
