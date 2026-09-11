# Oncogenicity Predictor Architecture Notes

## Current Direction

- Framework: FastAPI
- Deployment target: Render
- Public contract: single-variant `GET` and batch `POST`
- Response format: FHIR `Observation` for a single variant and a FHIR `Bundle` for batch results
- Scope: deterministic service implementation of the scoring approach described in https://pmc.ncbi.nlm.nih.gov/articles/PMC9081216/

## High-Level Flow

1. Validate incoming variant input.
2. Submit the variant to ClinGen for normalization.
3. Reject malformed or unnormalizable variants with a structured error response.
4. Use the normalized variant to collect annotations from VEP.
5. Run evidence pipelines.
6. Combine evidence into a final oncogenicity score.
7. Format the result and supporting evidence as FHIR.

## Evidence Pipelines

- Population data: primarily gnomAD
- Functional data: primarily MaveDB
- Predictive data: VEP and ClinVar
- Cancer hotspots
- Computational evidence: initially CADD

## Design Principles

- Normalize once, early, and keep a canonical internal variant representation.
- Keep source-specific client code separate from evidence interpretation.
- Keep scoring separate from HTTP and FHIR serialization.
- Treat partial evidence availability as a normal case rather than a fatal error.
- Preserve provenance for evidence and final scoring decisions.

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

## Small Starter Structure

```text
oncogenicity-predictor/
├── app/
│   ├── main.py
│   ├── api/
│   │   └── routes.py
│   ├── services/
│   │   ├── normalization/
│   │   │   └── variant_normalizer.py
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
│   └── models/
│       ├── requests.py
│       └── results.py
├── evaluation/
│   ├── datasets/
│   ├── runners/
│   ├── metrics/
│   └── reports/
├── tests/
├── docs/
│   └── architecture-notes.md
├── requirements.txt
├── render.yaml
└── README.md
```

## Near-Term Questions

1. What exact input format should the API require for submitted variants?
2. What are the minimum FHIR fields guaranteed in every response?
3. Which score output is authoritative for evaluation: numeric score, discrete class, or both?
4. Which pieces of `llm-oncogenicity` are worth migrating first?