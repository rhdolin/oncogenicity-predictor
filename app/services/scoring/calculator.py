from app.models.prediction import EvidenceResult, OncogenicityEvidence


def apply_evidence_interaction_rules(evidence: OncogenicityEvidence) -> OncogenicityEvidence:
    suppressed = evidence.model_copy(deep=True)
    evidence_results = (
        evidence.population,
        evidence.computational,
        evidence.hotspots,
        evidence.predictive,
        evidence.om1,
        evidence.op2,
        evidence.functional,
    )
    raw_codes = {
        result.evidenceCode
        for result in evidence_results
        if result.status != "not_available" and result.evidenceCode is not None
    }

    if "OS1" in raw_codes:
        suppressed.hotspots = _suppress_if_matches(
            suppressed.hotspots,
            target_code="OS3",
            reason="Suppressed because OS1 is applicable.",
        )

    if "OS1" in raw_codes or "OS3" in raw_codes:
        suppressed.om1 = _suppress_if_matches(
            suppressed.om1,
            target_code="OM1",
            reason="Suppressed because OS1 or OS3 is applicable.",
        )

    if "OS1" in raw_codes or "OS3" in raw_codes or "OM1" in raw_codes:
        suppressed.predictive = _suppress_if_matches(
            suppressed.predictive,
            target_code="OM4",
            reason="Suppressed because OS1, OS3, or OM1 is applicable.",
        )

    if "OM1" in raw_codes or "OM4" in raw_codes:
        suppressed.hotspots = _suppress_if_matches(
            suppressed.hotspots,
            target_code="OM3",
            reason="Suppressed because OM1 or OM4 is applicable.",
        )

    if "OVS1" in raw_codes:
        suppressed.predictive = _suppress_if_matches(
            suppressed.predictive,
            target_code="OM2",
            reason="Suppressed because OVS1 is applicable.",
        )

    return suppressed


def calculate_overall_score(evidence: OncogenicityEvidence) -> int:
    return (
        _counted_score(evidence.population)
        + _counted_score(evidence.computational)
        + _counted_score(evidence.hotspots)
        + _counted_score(evidence.predictive)
        + _counted_score(evidence.om1)
        + _counted_score(evidence.op2)
        + _counted_score(evidence.functional)
    )


def classify_overall_score(overall_score: int) -> str:
    if overall_score <= -7:
        return "Benign"
    if overall_score <= -1:
        return "Likely Benign"
    if overall_score <= 5:
        return "VUS"
    if overall_score <= 9:
        return "Likely Oncogenic"
    return "Oncogenic"


def _counted_score(evidence: EvidenceResult) -> int:
    if evidence.status != "applied":
        return 0
    return evidence.score


def _suppress_if_matches(
    evidence: EvidenceResult,
    target_code: str,
    reason: str,
) -> EvidenceResult:
    if evidence.evidenceCode != target_code or evidence.status != "applied":
        return evidence

    return evidence.model_copy(
        update={
            "status": "suppressed",
            "suppressionReason": reason,
            "dataAbsentReason": None,
        }
    )
