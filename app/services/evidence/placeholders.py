from app.models.prediction import EvidenceResult


def build_not_available_evidence(statement: str) -> EvidenceResult:
    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement=statement,
        status="not_available",
        source=None,
        matchedData=None,
        dataAbsentReason="unsupported",
    )
