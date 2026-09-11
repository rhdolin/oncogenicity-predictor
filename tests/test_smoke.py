from pathlib import Path
import sys

from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.main import app


client = TestClient(app)


def test_root_exposes_docs() -> None:
    response = client.get("/")

    assert response.status_code == 200
    assert response.json()["docs"] == "/docs"


def test_predict_single_returns_observation() -> None:
    response = client.get("/predict", params={"variant": "KRAS p.G12D"})

    assert response.status_code == 200
    body = response.json()
    assert body["resourceType"] == "Observation"
    assert body["subject"]["display"] == "KRAS p.G12D"


def test_predict_batch_returns_bundle() -> None:
    response = client.post("/predict", json={"variants": ["KRAS p.G12D", "BRAF p.V600E"]})

    assert response.status_code == 200
    body = response.json()
    assert body["resourceType"] == "Bundle"
    assert len(body["entry"]) == 2