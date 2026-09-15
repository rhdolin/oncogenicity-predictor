"""FastAPI application entrypoint for the oncogenicity predictor.

This file only owns app-level wiring such as metadata, router registration,
and lightweight service health endpoints. It does not contain variant
processing logic; that behavior lives in the route and service layers.
"""

from fastapi import FastAPI

from app.api.routes import router as api_router


app = FastAPI(
    title="Oncogenicity Predictor",
    description=(
        "Prototype API for oncogenicity prediction. The current implementation "
        "normalizes submitted variants through ClinGen, annotates them through "
        "Ensembl VEP, evaluates the currently implemented evidence pipelines, "
        "and exposes both an internal raw evidence summary endpoint plus FHIR "
        "Observation-style prediction responses. Submitted variants must "
        "currently be provided in HGVS format. Evidence and prediction endpoints "
        "also support an optional tumorType input for context-dependent rules. "
        "When VEP annotation fails, the prediction endpoints still return partial observations with "
        "component-level data absent reasons instead of failing the whole "
        "request."
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
    """Return a minimal service overview with a link to the interactive docs."""
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
    """Return a lightweight readiness response for uptime checks."""
    return {"status": "ok"}
