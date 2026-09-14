# Evidence And Scoring Specification

## Purpose

This document is the working specification for the oncogenicity evidence pipelines, the final score aggregation layer, and the current single-Observation prediction rendering.

- The architecture document describes where evidence and scoring fit in the system.
- This document defines the rule logic and response intent for each pipeline.
- For now, only the population pipeline is described in detail.
- The final client-facing prediction payload should not include `AnnotatedVariant` or `NormalizedVariant`.

## Scope

- Shared evidence result shape
- Population evidence rules
- Computational rules plus placeholders for hotspots, predictive, and functional pipelines
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

The exact enum values are now partially locked for the current implementation: `applied` and `not_available` are both in use, and the distinction between rule outcome and pipeline availability is preserved.

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

### Purpose

The current computational pipeline applies a narrow missense-only rule set using Ensembl VEP-exposed `CADD` and `FATHMM-XF` annotations.

### Current Upstream Inputs

The current annotation step populates:

- `computationalAnnotation.cadd.phred`
- `computationalAnnotation.cadd.raw`
- `computationalAnnotation.phyloP100wayVertebrate`
- `computationalAnnotation.fathmmXfCoding.prediction`
- `computationalAnnotation.fathmmXfCoding.score`
- `computationalAnnotation.fathmmXfCoding.rankscore`

Those values currently come from Ensembl REST VEP with:

- `CADD=true`
- `dbNSFP=ALL`

In the current REST response shape, the FATHMM-family fields exposed for this implementation are the `FATHMM-XF` dbNSFP keys with hyphenated names such as `fathmm-xf_coding_pred`.
Although those field names are the ones we read from the response, the current Ensembl REST service returned `invalid_field` when they were requested explicitly, so the implementation uses `dbNSFP=ALL` and then extracts the needed keys from the response.

### Current Intended Rule Basis

This rule set is explicitly scoped to missense variants.

- `OP1` is used as a positive supporting rule when `CADD PHRED >= 15`
- `SBP1` is used as a benign supporting rule when both of the following are true:
- `mostSevereConsequence == "missense_variant"`
- `CADD PHRED < 15`
- `FATHMM-XF` prediction is benign or neutral

Low `CADD` alone does not trigger `SBP1`. If `CADD` is low but the `FATHMM-XF` call is missing or not benign/neutral, the pipeline returns a score of `0` with no evidence code.

### Initial Rule Mapping

- If annotation failed or computational inputs are unavailable, return `not_available` with score `0`
- If `mostSevereConsequence` is not `missense_variant`, return score `0` with no evidence code
- If `CADD PHRED >= 15`, return `OP1` with score `1`
- Else if `CADD PHRED < 15` and `FATHMM-XF` is benign or neutral, return `SBP1` with score `-1`
- Else return score `0` with no evidence code

### Draft Decision Logic

```text
if annotation failed:
    not_available
elif most_severe_consequence != "missense_variant":
    score = 0
    evidenceCode = null
elif cadd_phred is missing:
    not_available
elif cadd_phred >= 15:
    score = 1
    evidenceCode = "OP1"
elif fathmm_xf_prediction in {"N", "neutral", "benign", "tolerated"}:
    score = -1
    evidenceCode = "SBP1"
else:
    score = 0
    evidenceCode = null
```

### Draft Evidence Statements

For `OP1`:

```json
{
  "score": 1,
  "evidenceCode": "OP1",
  "evidenceStatement": "CADD supports oncogenicity for this missense variant (PHRED 25.3).",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "mostSevereConsequence": "missense_variant",
    "caddPhred": 25.3,
    "caddRaw": 4.12,
    "phyloP100wayVertebrate": 7.89,
    "fathmmXfCodingPrediction": "D",
    "fathmmXfCodingScore": 0.88,
    "fathmmXfCodingRankscore": 0.91
  }
}
```

For `SBP1`:

```json
{
  "score": -1,
  "evidenceCode": "SBP1",
  "evidenceStatement": "Concordant computational predictors support a benign effect for this missense variant (CADD PHRED 10.4; FATHMM-XF N).",
  "status": "applied",
  "source": "vep",
  "matchedData": {
    "mostSevereConsequence": "missense_variant",
    "caddPhred": 10.4,
    "caddRaw": 0.42,
    "phyloP100wayVertebrate": 7.89,
    "fathmmXfCodingPrediction": "N",
    "fathmmXfCodingScore": 0.12,
    "fathmmXfCodingRankscore": 0.08
  }
}
```

For an evaluated missense variant with no current computational code trigger:

```json
{
  "score": 0,
  "evidenceCode": null,
  "evidenceStatement": "Computational evidence did not meet current scoring criteria.",
  "status": "applied",
  "source": "vep"
}
```

For a non-missense variant:

```json
{
  "score": 0,
  "evidenceCode": null,
  "evidenceStatement": "Computational missense rules were not applicable because the most severe consequence was synonymous_variant.",
  "status": "applied",
  "source": "vep"
}
```

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
- `Observation.extension`: custom extension carrying the originally submitted variant string in `valueString`
- `Observation.valueInteger`: overall numeric score
- `Observation.interpretation`: overall classification when available
- one `Observation.component` per evidence pipeline
- `component.code`: temporary code identifying the pipeline
- `component.valueInteger`: pipeline score when the pipeline is available
- `component.interpretation.coding.code`: pipeline evidence code such as `OP4`, `SBS1`, or `SBVS1`
- `component.interpretation.text`: short clinician-facing evidence statement
- `component.dataAbsentReason`: present instead of `component.valueInteger` when a pipeline is unavailable

Current serialization note:

- prediction endpoints omit `null` fields from the JSON response
- defaulted fields such as `resourceType="Observation"` and `status="final"` are still emitted

The final client-facing result should not embed `AnnotatedVariant` or `NormalizedVariant`.

Any normalization or annotation provenance that is needed by the client should be surfaced through the evidence output itself, primarily through fields such as:

- `source`
- `matchedData`
- `evidenceStatement`

If additional provenance is needed later, it should be added to the evidence result shape rather than by embedding internal intermediate objects into the final client-facing payload.

The exact aggregation rules are not yet locked in this document.

## Worked Examples

Placeholder. Add worked examples here drawn from current population and computational test fixtures, then expand them as additional pipelines land.