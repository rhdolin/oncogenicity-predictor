# Oncogenicity Predictor

This project predicts the oncogenicity of genomic variants and reports the result using the ClinGen/CGC/VICC oncogenicity framework described at https://pmc.ncbi.nlm.nih.gov/articles/PMC9081216/.

This repository contains a minimal FastAPI service that accepts HGVS variants, normalizes and annotates them, evaluates the currently implemented evidence pipelines, and returns an initial FHIR Observation-style oncogenicity prediction.

The live API can be explored and tested directly at https://oncogenicity-predictor.onrender.com/docs.

## Documentation

- [docs/evidence-and-scoring.md](docs/evidence-and-scoring.md): current evidence rule logic, data dependencies, and scoring semantics
- [docs/architecture-notes.md](docs/architecture-notes.md): current system architecture, flow, and implementation notes
- [docs/evaluation.md](docs/evaluation.md): evaluation workflow that compares the prediction output against a curated reference set.

## Endpoints

- `GET /`
- `GET /health`
- `GET /annotateVariant?variant=NM_004119.3%3Ac.2073T%3EG`
- `GET /summarizeEvidence?variant=NM_004119.3%3Ac.2073T%3EG[&tumorType=Breast%20Cancer]`
- `GET /predictOncogenicity?variant=NM_004119.3%3Ac.2073T%3EG[&tumorType=Breast%20Cancer]`
- `POST /predictOncogenicity`
- `GET /docs`

The current `GET /annotateVariant` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, annotates it through Ensembl VEP, and returns the internal `AnnotatedVariant` JSON shape.
The current `GET /summarizeEvidence` endpoint accepts a single variant in HGVS format and returns the raw evidence summary JSON before FHIR Observation mapping. It also accepts an optional `tumorType` query parameter used by context-dependent evidence rules. If no evidence lanes are evaluable, the summary returns a top-level `dataAbsentReason` plus a `predictionStatement`, and omits final score/classification.
The current `GET /predictOncogenicity` endpoint accepts a single variant in HGVS format, normalizes it through ClinGen, annotates it through Ensembl VEP, evaluates the currently implemented evidence pipelines, and returns a single FHIR Observation-style prediction object with final classification when available plus only score-contributing evidence components. If no evidence lanes are evaluable, the FHIR Observation instead uses a top-level `dataAbsentReason`. It also accepts an optional `tumorType` query parameter used by context-dependent evidence rules.
The current `POST /predictOncogenicity` endpoint accepts a list of variants in HGVS format and returns a FHIR `Bundle` containing one `Observation` resource per input variant. It also accepts an optional top-level `tumorType` field applied to every variant in the batch.
If VEP annotation fails for a variant, the prediction endpoints still return an in-band observation rather than failing the whole request, but the clinician-facing FHIR component list remains compact and may therefore be empty, with a top-level `dataAbsentReason` describing the unavailable overall prediction. The full audit surface lives in `GET /summarizeEvidence`.

## Local Run

```bash
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload
```

Then open `http://127.0.0.1:8000/docs`.

## Evaluation

The repository also includes a local evaluation workflow under `evaluation/`.

- `evaluation/variantLists/variantList.csv`: selected reference set used for evaluation
- `evaluation/variantLists/variantListMaster.csv`: larger reference-set source pool
- `evaluation/runEvaluation.py`: single local script that runs the predictor, flattens FHIR output into `oncogenicityPredictions.csv`, stores the exact `/summarizeEvidence`-style JSON in `evidenceSummary`, and writes the two metric CSVs

Run it from the repository root:

```bash
source .venv/bin/activate
python evaluation/runEvaluation.py
```

It overwrites these stable output targets:

- `evaluation/output/oncogenicityPredictions.csv`
- `evaluation/output/metrics-category-based.csv`
- `evaluation/output/metrics-score-based.csv`

`variantList.csv` is the canonical evaluation input. Fully specified rows use integer `points` values plus JSON-list `criteria` values, but category-only rows are also allowed when a source provides only the final classification. Those rows leave `points` and `criteria` blank, preserve any legacy free-text note in `comments`, still participate in category-based evaluation, and are excluded from score- and criteria-based calculations. `variantListMaster.csv` uses the same mixed-detail schema, with `referenceDetailLevel` distinguishing fully specified rows from category-only rows.

`oncogenicityPredictions.csv` keeps one row per variant from `variantList.csv` and includes the original reference-set columns translated into `reference*` fields, final FHIR-derived `predicted*` fields, a boolean `predictionUnavailable`, and the full internal `evidenceSummary` JSON for discrepancy review.

## Render Deployment

This project includes a minimal `render.yaml` blueprint.

Render will:

- install dependencies with `pip install -r requirements.txt`
- start the API with `uvicorn app.main:app --host 0.0.0.0 --port $PORT`
- check service health at `/health`

After deployment, the interactive API docs should be available at `/docs` on the Render URL, currently https://oncogenicity-predictor.onrender.com/docs.

## Current Scope

This is an incremental implementation.

- Current behavior: ClinGen-backed normalization, Ensembl VEP annotation, population, computational, hotspot, predictive, OM1, OP2, and ClinMAVE-backed functional evidence scoring, and single-Observation prediction output for the prediction endpoints
- Planned later behavior: additional evidence pipelines and further scoring refinements

The current rule set is deliberately narrower than mature manual-curation frameworks. Many disease-specific, gene-specific, and expert-panel-specific caveats and special cases are still documented limitations in v1 rather than automated logic.

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

The current single-variant prediction payload from `GET /predictOncogenicity` includes:

- top-level `Observation.issued` for the prediction timestamp
- top-level `Observation.extension` carrying the submitted variant HGVS string
- top-level `Observation.valueInteger` for the overall score when available
- top-level `Observation.interpretation` for the final classification when available
- top-level `Observation.dataAbsentReason` when the overall prediction is unavailable
- one `component` per score-contributing applied evidence lane
- pipeline `component.valueInteger` for included pipeline scores
- pipeline `component.interpretation` for evidence code plus short evidence statement

The current batch prediction payload from `POST /predictOncogenicity` includes:

- top-level `Bundle.resourceType = "Bundle"`
- top-level `Bundle.type = "collection"`
- top-level `Bundle.total` for the number of returned prediction resources
- one `entry.resource` per input variant, where each resource is the same FHIR `Observation` shape returned by `GET /predictOncogenicity`

Current functional evidence behavior includes:

- local ClinMAVE-backed lookup from `data/clinmave/variants.<GENE>.csv`
- lazy per-gene CSV loading with in-process caching on first use
- exact transcript-HGVS matching against `normalizedVariant.transcript_hgvs.mane_select_b38`
- `OS2` for oncogene plus gain-of-function and tumor suppressor gene plus loss-of-function
- `SBS2` score `-4` for functionally normal variants in resolved oncogene or tumor suppressor gene contexts
- score `0` with `applied` status for opposite-direction abnormal functional results
- score `0` with `applied` status when ClinMAVE contains conflicting functional classifications for the same exact matched variant
- `not_available` when ClinMAVE does not provide usable evidence for the queried gene, variant, or tumor-type context
- optional `tumorType` support on evidence and prediction endpoints for context-dependent genes such as `GATA3`

The current retained ClinMAVE gene panel is:

`AKT1`, `ALK`, `ARID1A`, `ATM`, `BAP1`, `BRAF`, `BRCA1`, `BRCA2`, `CDK4`, `CDK6`, `CHEK2`, `EGFR`, `ERBB2`, `GATA3`, `HRAS`, `JAK2`, `KRAS`, `MET`, `MYC`, `NF1`, `NRAS`, `PIK3CA`, `PTEN`, `SMAD4`, `TP53`, and `VHL`.

`data/clinmave/genes.txt` is the source-of-truth manifest for that retained panel.

The current computational pipeline is intentionally narrow:

- `OP1` is applied when `CADD PHRED >= 15`, including non-missense variants with usable CADD annotation
- `SBP1` is applied only for missense variants when `CADD PHRED < 15` and `FATHMM-XF` is concordantly benign or neutral
- `FATHMM-XF` values currently come from Ensembl REST VEP `dbNSFP` fields such as `fathmm-xf_coding_pred`

The bundled hotspot workbook at `data/hotspots_v2.xlsx` is sourced from Cancer Hotspots: https://www.cancerhotspots.org/#/home

The current hotspot pipeline loads that workbook into an in-memory cache on first use. It applies the legacy `OS3`, `OM3`, and `OP3` thresholds using exact gene plus protein-event matching, with SNVs requiring exact amino-acid substitution agreement and indels limited to direct matches supported by the workbook's native representation.

The current OM1 pipeline uses the curated local ClinGen-derived table at `data/om1_clingen_domains_seed.csv`. Runtime evaluation is intentionally limited to rows with `rowStatus=ready`, requires a MANE Select transcript consequence with a localized protein-altering event and resolvable residue position or span, and applies `OM1` with score `2` when that event overlaps a curated critical domain interval.

Coordinate conventions currently used by the normalizer include:

- `chrom` values like `chr13`, `chrX`, `chrY`, and `chrM`
- `chrom_num` values like `13`, `23` for X, `24` for Y, and `M` for mitochondrial variants

## Disclaimer

This repository is a rapid prototyping implementation designed to support experimentation.
It is a wholly automated algorithm based on curated knowledge sets that are evolving and generally capture only a subset of the knowledge present in the medical literature.
As a result, the algorithm can undercall oncogenicity scores, especially relative to pipelines that include manual curation.
It is not fit for actual clinical use and must not be used for patient care or clinical decision-making.