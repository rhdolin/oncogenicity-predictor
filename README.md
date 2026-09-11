# Oncogenicity Predictor

This repository currently contains a minimal FastAPI service for validating local startup, Render deployment, and the first normalization slice of the oncogenicity pipeline.

## Endpoints

- `GET /`
- `GET /health`
- `GET /predictOncogenicity?variant=NM_004119.3%3Ac.2073T%3EG`
- `POST /predictOncogenicity`
- `GET /docs`

The current `GET /predictOncogenicity` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, and returns the internal `NormalizedVariant` JSON shape.
The current `POST /predictOncogenicity` endpoint accepts a list of variants in HGVS format and returns a list of normalized variants in the same internal shape.

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

## Current Scope

This is an incremental implementation.

- Current behavior: ClinGen-backed normalization and internal normalized-variant output
- Planned later behavior: evidence pipelines, scoring, and final FHIR `Observation` or `Bundle` responses

At the moment, the normalization path requires submitted variants to be in HGVS format.

The current normalized payload can include:

- `geneSymbol`
- `geneNCBI_id`
- RefSeq-backed genomic, transcript, protein, and coordinate representations

Coordinate conventions currently used by the normalizer include:

- `chrom` values like `chr13`, `chrX`, `chrY`, and `chrM`
- `chrom_num` values like `13`, `23` for X, `24` for Y, and `M` for mitochondrial variants