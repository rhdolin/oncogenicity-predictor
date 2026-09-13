from fastapi import FastAPI

from app.api.routes import router as api_router


app = FastAPI(
    title="Oncogenicity Predictor",
    description=(
        "Prototype API for oncogenicity prediction. The current implementation "
        "normalizes submitted variants through ClinGen, annotates them through "
        "Ensembl VEP, and returns a single FHIR Observation-style prediction "
        "response for the prediction endpoints. Submitted variants must "
        "currently be provided in HGVS format. When VEP annotation fails, the "
        "prediction endpoints still return partial observations instead of "
        "failing the whole request."
    ),
    version="0.1.0",
)

app.include_router(api_router)


@app.get(
    "/",
    summary="Service overview",
    description="Returns a minimal status payload with a link to the interactive API docs.",
)
def root() -> dict:
    return {
        "service": "oncogenicity-predictor",
        "status": "ok",
        "docs": "/docs",
    }


@app.get(
    "/health",
    summary="Health check",
    description="Returns a lightweight readiness response for uptime checks and Render health probes.",
)
def health() -> dict:
    return {"status": "ok"}
