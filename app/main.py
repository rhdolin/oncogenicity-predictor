from fastapi import FastAPI
from pydantic import BaseModel, Field


class BatchRequest(BaseModel):
    variants: list[str] = Field(min_length=1)


def build_observation(variant: str) -> dict:
    return {
        "resourceType": "Observation",
        "status": "final",
        "code": {
            "text": "Predicted oncogenicity"
        },
        "valueCodeableConcept": {
            "text": "stub"
        },
        "subject": {
            "display": variant
        },
        "note": [
            {
                "text": "Stub response for deployment validation."
            }
        ]
    }


app = FastAPI(
    title="Oncogenicity Predictor",
    description="Stub service for Render deployment validation.",
    version="0.1.0",
)


@app.get("/")
def root() -> dict:
    return {
        "service": "oncogenicity-predictor",
        "status": "ok",
        "docs": "/docs",
    }


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.get("/predictOncogenicity")
def predict_single(variant: str) -> dict:
    return build_observation(variant)


@app.post("/predictOncogenicity")
def predict_batch(request: BatchRequest) -> dict:
    return {
        "resourceType": "Bundle",
        "type": "collection",
        "entry": [
            {"resource": build_observation(variant)}
            for variant in request.variants
        ],
    }