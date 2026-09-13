# Oncogenicity Predictor

This repository currently contains a minimal FastAPI service for validating local startup, Render deployment, and the first normalization-plus-annotation slice of the oncogenicity pipeline.

## Endpoints

- `GET /`
- `GET /health`
- `GET /annotateVariant?variant=NM_004119.3%3Ac.2073T%3EG`
- `GET /predictOncogenicity?variant=NM_004119.3%3Ac.2073T%3EG`
- `POST /predictOncogenicity`
- `GET /docs`

The current `GET /annotateVariant` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, annotates it through Ensembl VEP, and returns the internal `AnnotatedVariant` JSON shape.
The current `GET /predictOncogenicity` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, annotates it through Ensembl VEP, and returns the internal `AnnotatedVariant` JSON shape.
The current `POST /predictOncogenicity` endpoint accepts a list of variants in HGVS format and returns a list of annotated variants in the same internal shape.
If VEP annotation fails for a variant, the API still returns that variant's `AnnotatedVariant` with `annotationStatus="failed"`, `annotationError` populated, and annotation payload sections set to `null`.

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

- Current behavior: ClinGen-backed normalization, Ensembl VEP annotation, and internal annotated-variant output
- Planned later behavior: evidence pipelines, scoring, and final FHIR `Observation` or `Bundle` responses

At the moment, the normalization and annotation path requires submitted variants to be in HGVS format.

The current annotated payload includes:

- `normalizedVariant`, embedding the internal `NormalizedVariant`
- `annotationStatus` and `annotationError` for uniform success/failure reporting
- `basicAnnotation.mostSevereConsequence`
- RefSeq-only `basicAnnotation.transcriptConsequences`
- summarized `basicAnnotation.population`
- `computationalAnnotation.cadd`
- `computationalAnnotation.phyloP100wayVertebrate`

The embedded normalized payload can include:

- `geneSymbol`
- `geneNCBI_id`
- RefSeq-backed genomic, transcript, protein, and coordinate representations

Current annotation behavior includes:

- VEP query order: `genomic_hgvs.GRCh38`, then `transcript_hgvs.mane_select_b38`
- no additional VEP fallback after those two attempts in v1
- backfilling `transcript_hgvs.mane_select_b38` from VEP when ClinGen does not provide MANE Select but VEP does

Coordinate conventions currently used by the normalizer include:

- `chrom` values like `chr13`, `chrX`, `chrY`, and `chrM`
- `chrom_num` values like `13`, `23` for X, `24` for Y, and `M` for mitochondrial variants