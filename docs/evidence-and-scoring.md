# Evidence And Scoring Specification

## Purpose

This document is the working specification for the oncogenicity evidence pipelines and the final score aggregation layer.

- The architecture document describes where evidence and scoring fit in the system.
- This document defines the rule logic and response intent for each pipeline.
- For now, only the population pipeline is described in detail.
- The final client-facing prediction payload should not include `AnnotatedVariant` or `NormalizedVariant`.

## Scope

- Shared evidence result shape
- Population evidence rules
- Placeholders for computational, hotspots, predictive, and functional pipelines
- Placeholder for final score aggregation and classification mapping

## Shared Evidence Result Shape

Each pipeline is expected to produce a structured result shaped roughly like:

```json
{
  "score": 1,
  "evidenceCode": "OP4",
  "evidenceStatement": "Present at low frequency in gnomAD (≤1%; observed 0.20%).",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "maxSubpopulationAf": 0.002,
    "maxSubpopulationLabel": "gnomade_nfe",
    "maxOverallAf": 0.0015,
    "maxOverallLabel": "gnomadg"
  }
}
```

Working meanings for the additional fields:

- `status`: whether the rule was applied, evaluated but not triggered, not available, or failed
- `source`: the immediate source used by the pipeline, such as `vep`, `clinvar`, `cancerhotspots`, or `mavedb`
- `matchedData`: the structured values that drove the result so the rule is auditable without parsing free text

The exact enum values are not yet locked, but this document assumes the distinction between rule outcome and pipeline availability will be preserved.

## Population Pipeline

### Purpose

The population pipeline interprets gnomAD allele-frequency data and emits one evidence result for population-based oncogenicity support or benign support.

### Current Upstream Inputs

The current predictor implementation already collapses VEP population output into:

- `basicAnnotation.population.maxSubpopulationAf`
- `basicAnnotation.population.maxSubpopulationLabel`
- `basicAnnotation.population.maxOverallAf`
- `basicAnnotation.population.maxOverallLabel`

These values are derived from VEP co-located variant frequencies sourced from gnomAD.

More specifically, the current implementation does not query a standalone gnomAD release directly. It calls the current-assembly Ensembl REST VEP endpoint at `https://rest.ensembl.org/vep/human/hgvs`, which is the GRCh38 service in this repository, and it uses the gnomAD frequency fields that VEP returns inside `colocated_variants[].frequencies`. In v1 that means:

- gnomAD exomes overall and subpopulation fields: `gnomade` and `gnomade_*`
- gnomAD genomes overall and subpopulation fields: `gnomadg` and `gnomadg_*`

The current code takes the maximum relevant value across both of those VEP-exposed sources. It does not currently use a separate gnomAD structural-variant dataset.

As of the current Ensembl release documentation checked for this project, the human VEP data tables list `gnomAD exomes` as `v4.1` and `gnomAD genomes` as `v4.1`. That same published table is shown for both the main GRCh38 Ensembl release site and the GRCh37 Ensembl release site, so this repository should not document a simple `GRCh37 -> gnomAD v2.1.1` rule without stronger evidence.

The exact underlying gnomAD release is therefore still controlled by the Ensembl VEP service release we are calling, not by a separate configuration in this repository. This repository does not currently pin or configure a standalone gnomAD version independently of VEP.

For v1, the population pipeline should rely only on these collapsed summary fields rather than carrying forward the full raw allele-frequency map.

### Current Intended Rule Basis

The starting point is the existing deterministic population-frequency logic from the prior prototype implementation.

That logic:

1. Prefers the ALT allele bucket from VEP `allele_string` and `frequencies`.
2. Computes the maximum observed subpopulation allele frequency across `gnomade_*` and `gnomadg_*` population fields.
3. Separately computes the maximum observed overall allele frequency across `gnomade` and `gnomadg`.
4. Uses the maximum subpopulation allele frequency if available.
5. Falls back to the maximum overall allele frequency if no subpopulation value is available.
6. Falls back to `0.0` if VEP annotation succeeds but no usable population frequency is available.

In other words, the scoring rule is driven by a single effective population allele frequency chosen in this order:

1. `maxSubpopulationAf`
2. `maxOverallAf`
3. `0.0`

If VEP returns successfully but population-frequency fields are absent, the effective population allele frequency is treated as `0.0`.

### Initial Rule Mapping

The current intended rule mapping is:

- If effective population AF is greater than `0.05`, return `SBVS1` with score `-8`
- Else if effective population AF is greater than `0.01`, return `SBS1` with score `-4`
- Else return `OP4` with score `1`

This means the current draft behavior treats both of these as `OP4`:

- complete absence from controls
- low frequency less than or equal to `1%`

This policy is now locked for v1.

### Draft Decision Logic

Pseudo-logic for the population pipeline:

```text
effective_af = maxSubpopulationAf if present
    else maxOverallAf if present
    else 0.0

if effective_af > 0.05:
    score = -8
    evidenceCode = "SBVS1"
elif effective_af > 0.01:
    score = -4
    evidenceCode = "SBS1"
elif effective_af == 0.0:
    score = 1
    evidenceCode = "OP4"
else:
    score = 1
    evidenceCode = "OP4"
```

### Draft Evidence Statements

For `SBVS1`:

```json
{
  "score": -8,
  "evidenceCode": "SBVS1",
  "evidenceStatement": "Minor allele frequency is >5% in gnomAD (5.40%).",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "effectiveAf": 0.054,
    "effectiveAfSource": "maxSubpopulationAf",
    "maxSubpopulationAf": 0.054,
    "maxSubpopulationLabel": "gnomade_nfe",
    "maxOverallAf": 0.041,
    "maxOverallLabel": "gnomadg"
  }
}
```

For `SBS1`:

```json
{
  "score": -4,
  "evidenceCode": "SBS1",
  "evidenceStatement": "Minor allele frequency is >1% in gnomAD (1.40%).",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "effectiveAf": 0.014,
    "effectiveAfSource": "maxSubpopulationAf",
    "maxSubpopulationAf": 0.014,
    "maxSubpopulationLabel": "gnomadg_amr",
    "maxOverallAf": 0.009,
    "maxOverallLabel": "gnomade"
  }
}
```

For `OP4` when absent from controls:

```json
{
  "score": 1,
  "evidenceCode": "OP4",
  "evidenceStatement": "Absent from gnomAD controls.",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "effectiveAf": 0.0,
    "effectiveAfSource": "none",
    "maxSubpopulationAf": null,
    "maxSubpopulationLabel": null,
    "maxOverallAf": null,
    "maxOverallLabel": null
  }
}
```

For `OP4` when present at low frequency:

```json
{
  "score": 1,
  "evidenceCode": "OP4",
  "evidenceStatement": "Present at low frequency in gnomAD (≤1%; observed 0.20%).",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "effectiveAf": 0.002,
    "effectiveAfSource": "maxSubpopulationAf",
    "maxSubpopulationAf": 0.002,
    "maxSubpopulationLabel": "gnomade_nfe",
    "maxOverallAf": 0.0015,
    "maxOverallLabel": "gnomadg"
  }
}
```

### Handling Missing Annotation Data

If VEP annotation failed and population data is unavailable, the population pipeline should not pretend the rule was evaluated.

Current recommendation:

```json
{
  "score": 0,
  "evidenceCode": null,
  "evidenceStatement": "Population evidence could not be evaluated because annotation data was unavailable.",
  "status": "not_available",
  "source": "vep",
  "matchedData": null
}
```

This behavior is specifically for VEP annotation failure or lack of any VEP result. If VEP returns successfully but no population frequency data is present, the population pipeline should treat the effective allele frequency as `0.0` rather than `not_available`.

## Computational Pipeline

Placeholder. Expected initial inputs are `computationalAnnotation.cadd` and `computationalAnnotation.phyloP100wayVertebrate`.

## Hotspots Pipeline

Placeholder. Source selection and matching semantics are not yet locked.

## Predictive Pipeline

Placeholder. The split between consequence-based logic and ClinVar-based logic is not yet locked.

## Functional Pipeline

Placeholder. MaveDB is the expected first deterministic source, with any literature-based synthesis deferred until explicitly designed.

## Final Score Aggregation

The current API rendering for prediction endpoints is a single FHIR Observation-style object per variant.

Current rendering intent:

- `Observation.code`: temporary code for oncogenicity prediction
- `Observation.issued`: timestamp when the service generated the prediction
- `Observation.derivedFrom`: temporary reference carrying the originally submitted variant string
- `Observation.valueInteger`: overall numeric score
- `Observation.interpretation`: overall classification when available
- one `Observation.component` per evidence pipeline
- `component.code`: temporary code identifying the pipeline
- `component.valueInteger`: pipeline score when the pipeline is available
- `component.interpretation.coding.code`: pipeline evidence code such as `OP4`, `SBS1`, or `SBVS1`
- `component.interpretation.text`: short clinician-facing evidence statement
- `component.dataAbsentReason`: present instead of `component.valueInteger` when a pipeline is unavailable

The final client-facing result should not embed `AnnotatedVariant` or `NormalizedVariant`.

Any normalization or annotation provenance that is needed by the client should be surfaced through the evidence output itself, primarily through fields such as:

- `source`
- `matchedData`
- `evidenceStatement`

If additional provenance is needed later, it should be added to the evidence result shape rather than by embedding internal intermediate objects into the final client-facing payload.

The exact aggregation rules are not yet locked in this document.

## Worked Examples

Placeholder. Once the population pipeline is implemented, add examples here drawn from test fixtures.