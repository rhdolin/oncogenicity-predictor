# Oncogenicity Predictor Architecture Notes

## Current Direction

- Framework: FastAPI
- Deployment target: Render
- Public contract: single-variant `GET` endpoints and batch `POST`
- Current response format: `GET /annotateVariant` returns internal `AnnotatedVariant`; `GET /summarizeEvidence` returns the internal evidence summary; prediction endpoints return a single FHIR Observation-style prediction object
- Target later response format: keep the single-variant FHIR `Observation` shape and add a FHIR `Bundle` for batch results
- Target client-facing prediction payload: FHIR Observation content built from evidence summary plus final score/classification, without embedding internal `AnnotatedVariant` or `NormalizedVariant` objects
- Scope: deterministic service implementation of the scoring approach described in https://pmc.ncbi.nlm.nih.gov/articles/PMC9081216/

## High-Level Flow

1. Validate incoming HGVS variant input.
2. Submit the variant to ClinGen for normalization.
3. Reject malformed or unnormalizable variants with a structured error response instead of creating a `NormalizedVariant`.
4. Use the normalized variant to collect annotations from VEP.
5. Run evidence pipelines.
6. Combine evidence into a final oncogenicity score.
7. Format the result and supporting evidence as FHIR.

## Current Implementation Slice

- The deployed API has already been validated on Render.
- The current non-stub implementation slice is variant normalization, first-pass annotation, the population, computational, hotspot, predictive, OM1, OP2, and functional evidence pipelines, plus final score aggregation, interaction suppression, and compact FHIR Observation rendering.
- `GET /annotateVariant` returns the internal annotation-layer result, `GET /summarizeEvidence` returns the raw evidence summary before FHIR mapping, and `GET /predictOncogenicity` plus `POST /predictOncogenicity` return FHIR Observation-style prediction payloads.
- Submitted variants must currently be provided in HGVS format.
- The route layer calls orchestration entrypoints. The annotation-only flow delegates to the ClinGen-backed variant normalizer and then the VEP-backed variant annotator, while the prediction flow continues through evidence building, score aggregation, and FHIR Observation mapping.
- `canonical_b37` is currently populated only from a transcript allele that has a `genomeAlignments` entry for `GRCh37`, using the first primary `NM_` HGVS string from that transcript.
- `representative_transcript_hgvs` is currently populated from the best available NCBI RefSeq transcript in this order: `mane_select_b38`, then `canonical_b37`, then the first `NM_` transcript returned by ClinGen.
- `mane_select_b38` is populated from ClinGen when available; if ClinGen does not provide MANE Select, the annotation step can backfill it from VEP `hgvsc` on the MANE-marked RefSeq transcript row and mark `mane_select_b38_source` as `vep`.
- Coordinate normalization currently emits `chrM` for mitochondrial variants, and uses `23` for X plus `24` for Y in `chrom_num`.
- The VEP query strategy for v1 is intentionally narrow: try `genomic_hgvs.GRCh38` first, then `transcript_hgvs.mane_select_b38`, and stop there.

## Current AnnotatedVariant Slice

- `normalizedVariant` embeds the current internal normalization model without renaming its fields.
- `annotationStatus` is currently `complete` or `failed`.
- `annotationError` is present when VEP annotation fails and currently captures the source, message, and attempted query forms.
- `basicAnnotation.mostSevereConsequence` comes from VEP `most_severe_consequence`.
- `basicAnnotation.transcriptConsequences` keeps only RefSeq transcript rows whose `transcript_id` starts with `NM_`.
- Each retained transcript consequence currently includes `transcriptRefSeq`, `consequenceTerms`, `proteinStart`, `proteinEnd`, `aminoAcids`, `proteinHgvs`, `proteinEventType`, `rawProteinHgvs`, and `isManeSelect`.
- `basicAnnotation.population` collapses co-located allele frequencies into `maxSubpopulationAf`, `maxSubpopulationLabel`, `maxOverallAf`, and `maxOverallLabel`.
- `computationalAnnotation` is intentionally lean in v1 and currently includes `cadd`, `phyloP100wayVertebrate`, and `fathmmXfCoding`.
- On annotation failure, the API still returns `AnnotatedVariant` with `normalizedVariant` populated, `annotationStatus="failed"`, `annotationError` populated, and both annotation sections set to `null`.

## Evidence Pipelines

- Population data: Ensembl VEP co-located variant frequencies from gnomAD exomes (`gnomade*`) and gnomAD genomes (`gnomadg*`)
- Functional data: retained local ClinMAVE per-gene CSV exports
- OM1 domain data: curated local ClinGen-derived MANE-anchored CSV intervals
- Predictive data: VEP and ClinVar
- Cancer hotspots
- Computational evidence: broad `OP1` support from high `CADD`, plus missense-only `SBP1` benign concordance with `FATHMM-XF`

Functional evidence currently uses exact MANE transcript HGVS matching against ClinMAVE `Identifier` values after local normalization of the ClinMAVE identifier string. A lazy in-memory per-gene cache is used so the service reads one ClinMAVE CSV on first use rather than preloading the full retained panel at startup.

The current retained ClinMAVE panel is:

- `BRAF`
- `KRAS`
- `NRAS`
- `HRAS`
- `EGFR`
- `ERBB2`
- `ALK`
- `MET`
- `PIK3CA`
- `AKT1`
- `PTEN`
- `TP53`
- `NF1`
- `ARID1A`
- `SMAD4`
- `JAK2`
- `BRCA1`
- `BRCA2`
- `ATM`
- `CHEK2`
- `VHL`
- `BAP1`
- `CDK4`
- `CDK6`
- `GATA3`
- `MYC`

`data/clinmave/genes.txt` is the source-of-truth manifest for that retained local panel.

Prediction and evidence-summary endpoints now also accept an optional `tumorType` input for context-dependent evidence logic. Current concrete uses are `GATA3`, which can resolve to oncogene or tumor suppressor gene behavior depending on tumor type, and `OP2`, which uses a small curated rules table in `data/op2_rules.csv` with optional MANE protein HGVS matching.

OM1 uses `data/om1_clingen_domains_seed.csv` as a local ClinGen-derived runtime table. The current implementation only operationalizes rows marked `rowStatus=ready`, requires a MANE Select transcript consequence, and matches localized protein residue positions or spans against curated domain intervals.

Final scoring currently applies deterministic interaction suppression before summing scores. Suppressed evidence remains visible in the internal summary with `status="suppressed"` and a top-level `suppressionReason`, but only `applied` evidence contributes to the final `overallScore` and final classification.

The detailed rule specification for these pipelines and the final score layer lives in `docs/evidence-and-scoring.md`.

Client-facing prediction responses should surface any needed provenance through the evidence summary itself rather than by embedding internal normalization or annotation models.

## Design Principles

- Normalize once, early, and keep a canonical internal variant representation.
- Keep source-specific client code separate from evidence interpretation.
- Keep scoring separate from HTTP and from the dedicated FHIR serialization layer.
- Treat partial evidence availability as a normal case rather than a fatal error.
- Preserve provenance for evidence and final scoring decisions.
- Only instantiate `NormalizedVariant` on successful normalization.

## Evaluation Plan

- Evaluation does not need to be surfaced through the public API.
- Evaluation is implemented as a local Python script rather than an HTTP-driven workflow.
- The selected reference set lives at `evaluation/variantLists/variantList.csv`, retains the legacy `variant`, `gene`, `source`, `classification`, `points`, `criteria`, and `comments` field names, adds `referenceDetailLevel`, and supports both fully specified rows and category-only rows. When present, `points` are canonical integers and `criteria` are JSON-list values.
- The larger source pool lives at `evaluation/variantLists/variantListMaster.csv`, uses the same mixed-detail schema, and preserves category-only source rows by marking them with `referenceDetailLevel = category_only`, leaving `points` plus `criteria` blank, and carrying the original free-text criteria note in `comments`.
- The evaluation script builds the internal summary first, stores that exact `/summarizeEvidence`-style JSON as `evidenceSummary`, then derives the clinician-facing FHIR Observation and flattens that into one comparison row per input variant.
- The flattened predictions are written to `evaluation/output/oncogenicityPredictions.csv` with `reference*` and `predicted*` columns plus `predictionUnavailable`.
- Metrics are written to two separate CSVs: `evaluation/output/metrics-category-based.csv` for overall classification and score concordance, and `evaluation/output/metrics-score-based.csv` for exact criteria-set agreement and per-swimlane score agreement.
- Category-only rows still participate in classification concordance. Rows without reference scores are excluded from score concordance, and rows without reference criteria are excluded from criteria and swimlane concordance.
- Rows with `predictionUnavailable = true` are excluded from concordance metrics and counted separately.

## Reuse Plan From `llm-oncogenicity`

- Likely reusable: batch processing, evaluator logic, data-source adapters, parsing helpers, fixtures
- Likely to replace: LLM and RAG orchestration, older normalization flow, outdated MaveDB scoring assumptions, non-FHIR output surfaces
- Migration strategy: treat the old repository as a source of reusable modules rather than as the base architecture for the new service

## Current Evaluation Structure

```text
oncogenicity-predictor/
├── evaluation/
│   ├── variantLists/
│   │   ├── variantList.csv
│   │   └── variantListMaster.csv
│   ├── output/
│   │   ├── oncogenicityPredictions.csv
│   │   ├── metrics-category-based.csv
│   │   └── metrics-score-based.csv
│   └── runEvaluation.py
```

## Near-Term Questions

1. How strict should HGVS validation become before calling ClinGen?
2. Which additional FHIR fields beyond the current Observation skeleton should be locked before more pipelines land?
3. Which score output is authoritative for evaluation: numeric score, discrete class, or both?
4. Which pieces of `llm-oncogenicity` are worth migrating first?

## NormalizedVariant Shape

The current internal normalization target is intentionally permissive.

- Required: `submitted_variant`
- Required: `normalization.source`
- Required: `normalization.queried_variant`
- Optional: `identifiers.caid`
- Optional: `geneSymbol`
- Optional: `geneNCBI_id`
- Optional: all genomic, transcript, protein, and coordinate representations
- Optional: `protein.civic_profile_name`, which is currently computed locally from `geneSymbol` plus the derived protein short name
- Optional: `transcript_hgvs.mane_select_b38_source`, currently `clingen` or `vep` when `mane_select_b38` is populated

Notes:

- `uniprot_id` is intentionally deferred and is not currently part of the normalization model.
- The current implementation assumes successful normalization returns a `NormalizedVariant`; malformed or unnormalizable input returns an error instead.
- Normalized genomic, transcript, and protein fields are populated only from NCBI RefSeq accessions: `NC_`, `NM_`, and `NP_`.
- `representative_transcript_hgvs` is available as a practical fallback when MANE and canonical transcript fields are absent.

## Current Public Response Shape

- `GET /annotateVariant` currently returns `AnnotatedVariant`.
- `GET /annotateVariant` is the explicit annotation-oriented single-variant endpoint.
- `GET /summarizeEvidence` currently returns `OncogenicityPredictionSummary` and is intended as an internal/debug endpoint.
- Prediction endpoints currently return a single FHIR Observation-style object with `issued`, a custom extension carrying the submitted variant HGVS string, an overall score, a final classification, and only the score-contributing evidence components.
- Batch prediction requests currently return an `observations` list of those prediction objects.
- Prediction success/failure is currently expressed through the internal summary surface rather than by embedding the internal annotation result into the clinician-facing FHIR output.
- FHIR rendering currently lives in `app/services/fhir/observation_builder.py` and is intentionally lightweight rather than profile-complete.

Deferred manuscript caveats worth future implementation are currently documented in the evidence-and-scoring spec rather than automated. The highest-value deferred items remain `OVS1` splice and 3' end nuance, splicing-aware suppression of protein-level criteria, hotspot caution for truncating-driven hotspots, functional evidence downgrading, and hereditary predisposition population-threshold overrides.

The object is only created on successful normalization. Failures are handled as errors rather than partial `NormalizedVariant` instances.