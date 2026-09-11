# Oncogenicity Predictor

This repository currently contains a minimal FastAPI service for validating local startup and Render deployment.

## Endpoints

- `GET /`
- `GET /health`
- `GET /predict?variant=KRAS%20p.G12D`
- `POST /predict`
- `GET /docs`

The current `GET /predict` endpoint returns a stub FHIR `Observation`.
The current `POST /predict` endpoint returns a stub FHIR `Bundle`.

## Local Run

```bash
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload
```

Then open `http://127.0.0.1:8000/docs`.

## Render Deployment

This project includes a minimal `render.yaml` blueprint.

Render will:

- install dependencies with `pip install -r requirements.txt`
- start the API with `uvicorn app.main:app --host 0.0.0.0 --port $PORT`
- check service health at `/health`

After deployment, the interactive API docs should be available at `/docs` on the Render URL.