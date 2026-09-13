from app.models.prediction import OncogenicityEvidence


def calculate_overall_score(evidence: OncogenicityEvidence) -> int:
    return (
        evidence.population.score
        + evidence.computational.score
        + evidence.hotspots.score
        + evidence.predictive.score
        + evidence.functional.score
    )
