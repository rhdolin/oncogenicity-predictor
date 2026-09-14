# Oncogenicity Predictor

This repository currently contains a minimal FastAPI service for validating local startup, Render deployment, normalization, annotation, the population, computational, and hotspot evidence pipelines, and an initial FHIR Observation-style prediction output.

## Endpoints

- `GET /`
- `GET /health`
- `GET /annotateVariant?variant=NM_004119.3%3Ac.2073T%3EG`
- `GET /summarizeEvidence?variant=NM_004119.3%3Ac.2073T%3EG`
- `GET /predictOncogenicity?variant=NM_004119.3%3Ac.2073T%3EG`
- `POST /predictOncogenicity`
- `GET /docs`

The current `GET /annotateVariant` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, annotates it through Ensembl VEP, and returns the internal `AnnotatedVariant` JSON shape.
The current `GET /summarizeEvidence` endpoint is an internal/debug endpoint that accepts a single variant in HGVS format and returns the raw evidence summary JSON before FHIR Observation mapping.
The current `GET /predictOncogenicity` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, annotates it through Ensembl VEP, evaluates the currently implemented evidence pipelines, and returns a single FHIR Observation-style prediction object.
The current `POST /predictOncogenicity` endpoint accepts a list of variants in HGVS format and returns an `observations` list in that same shape.
If VEP annotation fails for a variant, the prediction endpoints still return a partial observation with `component.dataAbsentReason` set for the unavailable pipeline rather than failing the whole request.

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

- Current behavior: ClinGen-backed normalization, Ensembl VEP annotation, population, computational, and hotspot evidence scoring, and single-Observation prediction output for the prediction endpoints
- Planned later behavior: additional evidence pipelines, richer scoring/classification, and batch FHIR `Bundle` responses

At the moment, the normalization and annotation path requires submitted variants to be in HGVS format.

The current annotation payload from `GET /annotateVariant` includes:

- `normalizedVariant`, embedding the internal `NormalizedVariant`
- `annotationStatus` and `annotationError` for uniform success/failure reporting
- `basicAnnotation.mostSevereConsequence`
- RefSeq-only `basicAnnotation.transcriptConsequences`
- summarized `basicAnnotation.population`
- `computationalAnnotation.cadd`
- `computationalAnnotation.phyloP100wayVertebrate`
- `computationalAnnotation.fathmmXfCoding`

The embedded normalized payload can include:

- `geneSymbol`
- `geneNCBI_id`
- RefSeq-backed genomic, transcript, protein, and coordinate representations

Current annotation behavior includes:

- VEP query order: `genomic_hgvs.GRCh38`, then `transcript_hgvs.mane_select_b38`
- no additional VEP fallback after those two attempts in v1
- backfilling `transcript_hgvs.mane_select_b38` from VEP when ClinGen does not provide MANE Select but VEP does

The current prediction payload from `GET /predictOncogenicity` and `POST /predictOncogenicity` includes:

- top-level `Observation.issued` for the prediction timestamp
- top-level `Observation.extension` carrying the submitted variant HGVS string
- top-level `Observation.valueInteger` for the overall score
- top-level `Observation.interpretation`, currently omitted because classification is not yet implemented
- one `component` per evidence pipeline
- pipeline `component.valueInteger` for available pipeline scores
- pipeline `component.interpretation` for evidence code plus short evidence statement
- pipeline `component.dataAbsentReason` for unavailable pipelines

The current computational pipeline is intentionally narrow:

- `OP1` is applied when `CADD PHRED >= 15`, including non-missense variants with usable CADD annotation
- `SBP1` is applied only for missense variants when `CADD PHRED < 15` and `FATHMM-XF` is concordantly benign or neutral
- `FATHMM-XF` values currently come from Ensembl REST VEP `dbNSFP` fields such as `fathmm-xf_coding_pred`

The bundled hotspot workbook at `data/hotspots_v2.xlsx` is sourced from Cancer Hotspots: https://www.cancerhotspots.org/#/home

The current hotspot pipeline loads that workbook into an in-memory cache on first use. It applies the legacy `OS3`, `OM3`, and `OP3` thresholds using exact gene plus protein-event matching, with SNVs requiring exact amino-acid substitution agreement and indels limited to direct matches supported by the workbook's native representation.

Coordinate conventions currently used by the normalizer include:

- `chrom` values like `chr13`, `chrX`, `chrY`, and `chrM`
- `chrom_num` values like `13`, `23` for X, `24` for Y, and `M` for mitochondrial variants