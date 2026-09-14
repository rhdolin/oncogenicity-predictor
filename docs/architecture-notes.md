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
- The current non-stub implementation slice is variant normalization, first-pass annotation, the population and computational evidence pipelines, and initial FHIR Observation rendering.
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
- Each retained transcript consequence currently includes `transcriptRefSeq`, `consequenceTerms`, `proteinStart`, `proteinEnd`, `aminoAcids`, and `isManeSelect`.
- `basicAnnotation.population` collapses co-located allele frequencies into `maxSubpopulationAf`, `maxSubpopulationLabel`, `maxOverallAf`, and `maxOverallLabel`.
- `computationalAnnotation` is intentionally lean in v1 and currently includes `cadd`, `phyloP100wayVertebrate`, and `fathmmXfCoding`.
- On annotation failure, the API still returns `AnnotatedVariant` with `normalizedVariant` populated, `annotationStatus="failed"`, `annotationError` populated, and both annotation sections set to `null`.

## Evidence Pipelines

- Population data: Ensembl VEP co-located variant frequencies from gnomAD exomes (`gnomade*`) and gnomAD genomes (`gnomadg*`)
- Functional data: primarily MaveDB
- Predictive data: VEP and ClinVar
- Cancer hotspots
- Computational evidence: broad `OP1` support from high `CADD`, plus missense-only `SBP1` benign concordance with `FATHMM-XF`

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
- A practical first approach is to use the batch `POST` response as the prediction artifact.
- The returned FHIR `Bundle` should be flattened into a simple comparison table.
- The flattened predictions can then be compared against a gold-standard CSV.
- Metrics should emphasize concordance and classification accuracy, with room for stratified analysis later.

## Reuse Plan From `llm-oncogenicity`

- Likely reusable: batch processing, evaluator logic, data-source adapters, parsing helpers, fixtures
- Likely to replace: LLM and RAG orchestration, older normalization flow, outdated MaveDB scoring assumptions, non-FHIR output surfaces
- Migration strategy: treat the old repository as a source of reusable modules rather than as the base architecture for the new service

## Illustrative Future Structure

```text
oncogenicity-predictor/
├── app/
│   ├── main.py
│   ├── api/
│   │   └── routes.py
│   ├── models/
│   │   ├── annotated_variant.py
│   │   ├── normalized_variant.py
│   │   └── requests.py
│   ├── services/
│   │   ├── normalization/
│   │   │   └── variant_normalizer.py
│   │   ├── orchestration/
│   │   │   └── single_variant_pipeline.py
│   │   ├── annotation/
│   │   │   └── variant_annotator.py
│   │   ├── evidence/
│   │   │   ├── population.py
│   │   │   ├── functional.py
│   │   │   ├── predictive.py
│   │   │   ├── hotspots.py
│   │   │   └── computational.py
│   │   ├── scoring/
│   │   │   └── calculator.py
│   │   └── fhir/
│   │       └── observation_builder.py
├── evaluation/
│   ├── datasets/
│   ├── runners/
│   ├── metrics/
│   └── reports/
├── tests/
├── docs/
│   └── architecture-notes.md
│   └── evidence-and-scoring.md
├── requirements.txt
├── render.yaml
└── README.md
```

This section is aspirational rather than a verbatim snapshot of the current repo layout.

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
- Prediction endpoints currently return a single FHIR Observation-style object with `issued`, a custom extension carrying the submitted variant HGVS string, an overall score, and one component per evidence pipeline.
- Batch prediction requests currently return an `observations` list of those prediction objects.
- Prediction success/failure is currently expressed through component-level values versus `dataAbsentReason`, rather than by embedding the internal annotation result.
- FHIR rendering currently lives in `app/services/fhir/observation_builder.py` and is intentionally lightweight rather than profile-complete.

The object is only created on successful normalization. Failures are handled as errors rather than partial `NormalizedVariant` instances.