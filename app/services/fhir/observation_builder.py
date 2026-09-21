"""Map the internal prediction summary into a lightweight FHIR Observation.

This module is intentionally only a projection layer. It does not evaluate any
evidence itself; instead it turns the already-scored internal summary into the
clinician-facing Observation shape used by the prediction endpoints.
"""

from datetime import datetime, timezone

from app.models.prediction import (
    Annotation,
    BundleEntry,
    CodeableConcept,
    Coding,
    EvidenceResult,
    Extension,
    ObservationComponent,
    OncogenicityObservation,
    OncogenicityObservationBundle,
    OncogenicityPredictionSummary,
)


TEMP_CODE_SYSTEM = "https://oncogenicity-predictor.example/fhir/CodeSystem/temp-codes"
DATA_ABSENT_REASON_SYSTEM = "http://terminology.hl7.org/CodeSystem/data-absent-reason"
SUBMITTED_VARIANT_EXTENSION_URL = (
    "https://oncogenicity-predictor.example/fhir/StructureDefinition/submitted-variant"
)
PROTOTYPE_DISCLAIMER = (
    "This oncogenicity predictor is a rapid prototyping implementation intended to support experimentation. "
    "It is a wholly automated algorithm based on curated knowledge sets that are evolving and generally capture only a subset of the knowledge present in the medical literature. "
    "As a result, the algorithm can undercall oncogenicity scores, especially relative to pipelines that include manual curation. "
    "This algorithm is not fit for actual clinical use and must not be used "
    "for patient care or clinical decision-making."
)
OBSERVATION_METHOD_TEXT = (
    "OncogenicityPredictor v1 (https://github.com/rhdolin/oncogenicity-predictor)"
)


def _build_concept(
    code: str | None = None,
    display: str | None = None,
    text: str | None = None,
) -> CodeableConcept:
    """Build a small CodeableConcept using the temporary code system used in this prototype."""
    coding = []
    if code is not None or display is not None:
        coding.append(Coding(system=TEMP_CODE_SYSTEM, code=code, display=display))
    return CodeableConcept(coding=coding, text=text)


def _build_pipeline_component(
    pipeline_name: str,
    evidence: EvidenceResult,
) -> ObservationComponent:
    """Map one internal evidence result into a single Observation component."""
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


def _should_include_pipeline_component(evidence: EvidenceResult) -> bool:
    return evidence.status == "applied" and evidence.score != 0


def _build_scored_components(
    summary: OncogenicityPredictionSummary,
) -> list[ObservationComponent]:
    components: list[ObservationComponent] = []
    for pipeline_name, evidence in (
        ("population", summary.oncogenicityEvidence.population),
        ("computational", summary.oncogenicityEvidence.computational),
        ("hotspots", summary.oncogenicityEvidence.hotspots),
        ("predictive", summary.oncogenicityEvidence.predictive),
        ("om1", summary.oncogenicityEvidence.om1),
        ("op2", summary.oncogenicityEvidence.op2),
        ("functional", summary.oncogenicityEvidence.functional),
    ):
        if _should_include_pipeline_component(evidence):
            components.append(_build_pipeline_component(pipeline_name, evidence))
    return components


def build_oncogenicity_observation(
    summary: OncogenicityPredictionSummary,
    submitted_variant: str,
) -> OncogenicityObservation:
    """Project the internal summary into the current single-Observation FHIR response."""
    interpretation: list[CodeableConcept] = []
    if summary.overallClassification is not None:
        interpretation = [
            _build_concept(
                code=summary.overallClassification,
                display=summary.overallClassification,
                text=summary.overallClassification,
            )
        ]

    observation = OncogenicityObservation(
        issued=datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
            "+00:00",
            "Z",
        ),
        code=_build_concept(
            code="oncogenicity-prediction",
            display="Oncogenicity prediction",
            text="Oncogenicity prediction",
        ),
        method=_build_concept(text=OBSERVATION_METHOD_TEXT),
        valueInteger=summary.overallScore,
        interpretation=interpretation,
        extension=[
            Extension(
                url=SUBMITTED_VARIANT_EXTENSION_URL,
                valueString=submitted_variant,
            )
        ],
        note=[Annotation(text=PROTOTYPE_DISCLAIMER)],
        component=_build_scored_components(summary),
    )

    if summary.dataAbsentReason is not None:
        observation.dataAbsentReason = CodeableConcept(
            coding=[
                Coding(
                    system=DATA_ABSENT_REASON_SYSTEM,
                    code=summary.dataAbsentReason,
                    display=summary.dataAbsentReason,
                )
            ],
            text=summary.predictionStatement,
        )

    return observation


def build_oncogenicity_observation_bundle(
    observations: list[OncogenicityObservation],
) -> OncogenicityObservationBundle:
    return OncogenicityObservationBundle(
        total=len(observations),
        entry=[BundleEntry(resource=observation) for observation in observations],
    )
