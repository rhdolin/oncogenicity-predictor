# Evidence And Scoring Specification

## Purpose

This document is the working specification for the oncogenicity evidence pipelines, the current overall score aggregation layer, and the current single-Observation prediction rendering.

- The architecture document describes where evidence and scoring fit in the system.
- This document defines the rule logic, data dependencies, and response semantics for each evidence pipeline.
- The final client-facing prediction payload should not include `AnnotatedVariant` or `NormalizedVariant`.

## Current State Summary

The repository currently implements these evidence pipelines:

- population
- computational
- hotspots
- predictive
- om1
- op2
- functional

Current overall scoring applies a small deterministic interaction-resolution step, sums only score-contributing evidence, maps the adjusted total into a final classification, and returns one FHIR Observation-style result per queried variant.

## Scope

- shared evidence result shape and semantics
- rule order and availability behavior for each implemented pipeline
- current data sources and matching strategies
- current implemented pipelines plus remaining planned v1 additions
- current overall score aggregation and prediction rendering

## Shared Evidence Result Shape

Each pipeline emits one `EvidenceResult`.

Current model shape:

```json
{
  "score": 1,
  "evidenceCode": "OP4",
  "evidenceStatement": "Present at low frequency in gnomAD (<=1%; observed 0.20%).",
  "status": "applied",
  "suppressionReason": null,
  "source": "vep",
  "matchedData": {
    "maxSubpopulationAf": 0.002,
    "maxSubpopulationLabel": "gnomade_nfe",
    "maxOverallAf": 0.0015,
    "maxOverallLabel": "gnomadg"
  },
  "dataAbsentReason": null
}
```

Current field meanings:

- `score`: signed integer contribution from the pipeline
- `evidenceCode`: the applied evidence code when one exists, otherwise `null`
- `evidenceStatement`: human-readable explanation of the outcome
- `status`: currently constrained to `applied`, `suppressed`, or `not_available`
- `suppressionReason`: top-level explanation when a matched criterion is excluded from final scoring
- `source`: immediate source or pipeline label used by the rule
- `matchedData`: structured audit trail for the decision when available
- `dataAbsentReason`: reason token used when the pipeline is not available

Current status semantics are intentionally narrow:

- `applied` means the pipeline had enough information to evaluate its current rule set and its score contributes to the final sum, including score `0` outcomes where no rule fired or where matching data supported a neutral outcome
- `suppressed` means the pipeline matched a criterion with its original score and code preserved, but that criterion is excluded from final scoring by a deterministic interaction rule
- `not_available` means the pipeline could not be defensibly evaluated because required annotation, context, or local source data was missing, unreadable, or unsupported

Each evidence pipeline emits at most one result. When a pipeline contains multiple criteria, they are evaluated in a defined priority order and the first matching criterion wins. If no criterion matches but the pipeline had enough information to evaluate its current rule set, the pipeline returns a single `applied` result with score `0`.

At the overall summary level, the service can also report prediction unavailability. When no evidence lanes are evaluable, the summary returns:

- `overallScore = null`
- `overallClassification = null`
- top-level `dataAbsentReason`
- top-level `predictionStatement`

## Population Pipeline

### Purpose

The population pipeline interprets gnomAD allele-frequency data and emits one evidence result for population-based oncogenicity support or benign support.

### Current Upstream Inputs

The current predictor implementation collapses VEP population output into:

- `basicAnnotation.population.maxSubpopulationAf`
- `basicAnnotation.population.maxSubpopulationLabel`
- `basicAnnotation.population.maxOverallAf`
- `basicAnnotation.population.maxOverallLabel`

These values are derived from Ensembl REST VEP `colocated_variants[].frequencies` fields.

More specifically, the repository does not query a standalone gnomAD release directly. It calls the current-assembly Ensembl REST VEP endpoint and uses the gnomAD frequency fields that VEP returns. In v1 that means:

- gnomAD exomes overall and subpopulation fields: `gnomade` and `gnomade_*`
- gnomAD genomes overall and subpopulation fields: `gnomadg` and `gnomadg_*`

The exact underlying gnomAD release is therefore controlled by the Ensembl VEP service release rather than by a separate repository-level configuration.

### Current Decision Basis

The current logic:

1. Prefers the ALT allele bucket from VEP `allele_string` and `frequencies`.
2. Computes the maximum observed subpopulation allele frequency across `gnomade_*` and `gnomadg_*` population fields.
3. Separately computes the maximum observed overall allele frequency across `gnomade` and `gnomadg`.
4. Uses the maximum subpopulation allele frequency if available.
5. Falls back to the maximum overall allele frequency if no subpopulation value is available.
6. Falls back to `0.0` if VEP annotation succeeds but no usable population frequency is available.

The scoring rule is driven by one effective allele frequency chosen in this order:

1. `maxSubpopulationAf`
2. `maxOverallAf`
3. `0.0`

### Current Rule Mapping

- effective AF `> 0.05` -> `SBVS1`, score `-8`
- else effective AF `> 0.01` -> `SBS1`, score `-4`
- else -> `OP4`, score `1`

This means the current v1 policy treats both of these as `OP4`:

- complete absence from controls
- low frequency less than or equal to `1%`

### Current Availability Semantics

- `not_available` when annotation failed and population data is unavailable
- `applied` when annotation succeeded, even if no population frequency fields were present, because that case is interpreted as effective AF `0.0`

## Computational Pipeline

### Purpose

The computational pipeline applies a narrow rule set using Ensembl VEP-exposed `CADD` and `FATHMM-XF` annotations.

### Current Upstream Inputs

The current annotation step populates:

- `computationalAnnotation.cadd.phred`
- `computationalAnnotation.cadd.raw`
- `computationalAnnotation.phyloP100wayVertebrate`
- `computationalAnnotation.fathmmXfCoding.prediction`
- `computationalAnnotation.fathmmXfCoding.score`
- `computationalAnnotation.fathmmXfCoding.rankscore`
- `basicAnnotation.mostSevereConsequence`

These values currently come from Ensembl REST VEP with:

- `CADD=true`
- `dbNSFP=ALL`

The implementation extracts the needed `FATHMM-XF` fields from the returned `dbNSFP` payload because explicit field requests proved unreliable for this repository.

### Current Rule Order

1. `OP1`
2. `SBP1`
3. else score `0` with `applied`

### Current Rule Mapping

- `OP1` when `CADD PHRED >= 15`, including non-missense variants
- `SBP1` when all of the following are true:
- `mostSevereConsequence == "missense_variant"`
- `CADD PHRED < 15`
- `FATHMM-XF` prediction is benign or neutral

Low `CADD` alone does not trigger `SBP1`. Non-missense variants are not eligible for the current benign computational rule.

### Current Availability Semantics

- `not_available` when annotation failed
- `not_available` when `CADD` is unavailable
- `applied` score `0` when computational data was sufficient to evaluate the implemented rule set but no code fired

## Hotspots Pipeline

### Purpose

The hotspots pipeline evaluates the bundled Cancer Hotspots workbook against transcript-level protein consequences.

### Current Data Source

- local workbook: `data/hotspots_v2.xlsx`
- source provenance: Cancer Hotspots

### Current Data Access Strategy

- the workbook is parsed into an in-memory index on first use
- the service does not reopen the workbook for every variant
- unreadable or missing workbook data is treated as pipeline unavailability rather than as a score `0` evaluation

### Current Match Policy

- matching is gene-centric because the workbook is not transcript-indexed
- transcript consequences are evaluated in order, and the first defensible hotspot match wins
- SNVs require exact match on gene, amino-acid position, reference amino acid, and alternate amino acid
- indels are supported only when the transcript consequence can be compared directly to the workbook's native event representation
- there is no separate indel normalization layer beyond the current direct event matching

### Current Rule Mapping

- `OS3` when hotspot `mutation_count >= 50` and exact protein-event `variant_count >= 10`, score `4`
- `OM3` when exact protein-event `variant_count >= 10` but `OS3` does not apply, score `2`
- `OP3` when exact protein-event `variant_count` is between `1` and `9`, score `1`
- else score `0`

### Current Availability Semantics

- `not_available` when annotation failed
- `not_available` when the normalized gene symbol is unavailable
- `not_available` when the local hotspot workbook is missing or unreadable
- `applied` score `0` when no transcript consequence yields a defensible hotspot match

## Predictive Pipeline

### Purpose

The predictive pipeline combines deterministic consequence-based logic with ClinVar-backed protein-change matching.

### Current Upstream Inputs

The current implementation depends on:

- `normalizedVariant.geneSymbol`
- `normalizedVariant.protein.short_name`
- `normalizedVariant.protein.hgvs_3letter`
- `basicAnnotation.mostSevereConsequence`
- `computationalAnnotation.phyloP100wayVertebrate`
- optional `tumorType`
- local gene-role metadata in `data/_Dict_Gene.csv`
- ClinVar E-utilities search and summary lookups

### Current Role Resolution Policy

Most genes resolve directly from `data/_Dict_Gene.csv` as oncogene, tumor suppressor gene, both, or neither.

`GATA3` is currently the main dual-role case. Predictive interpretation uses optional `tumorType` to resolve whether the queried context should behave as oncogene or tumor suppressor gene. If a dual-role gene cannot be resolved for the supplied tumor type, role-dependent predictive rules are treated as unavailable rather than evaluated as neutral.

### Current Rule Order

1. `OVS1`
2. `OS1`
3. `OM2`
4. `SBP2`
5. `OM4`
6. else score `0` with `applied`

### Current Rule Mapping

- `OVS1`, score `8`: null variant in a resolved tumor suppressor gene context
- `OS1`, score `4`: same amino acid change as a previously established somatic oncogenic ClinVar variant
- `OM2`, score `2`: in-frame insertion or deletion in a resolved oncogene or tumor suppressor gene context, or `stop_lost` in a resolved tumor suppressor gene context
- `SBP2`, score `-1`: synonymous variant with `phyloP100wayVertebrate < 2.0`
- `OM4`, score `2`: missense variant at an amino-acid residue where a different missense variant is established as somatic oncogenic in ClinVar

### Current ClinVar Match Policy

- ClinVar ESearch uses fielded queries built from gene plus protein tokens
- Entrez `[varnam]` search is not treated as exact matching
- returned ClinVar summaries are filtered locally against exact protein aliases before `OS1` or `OM4` is applied
- ClinVar lookup failure is treated as `not_available` for the affected rule path rather than as a neutral score `0`

### Current Availability Semantics

- `not_available` when annotation failed
- `not_available` when required local fields for a rule path are unavailable, such as missing gene symbol, missing residue token, or missing `phyloP`
- `not_available` when a dual-role gene requires tumor type but tumor type is missing or unresolved for the current mapping
- `not_available` when the relevant ClinVar lookup fails
- `applied` score `0` when the pipeline had enough information to evaluate the current predictive rule set but no rule fired

## Functional Pipeline

### Purpose

The functional pipeline uses locally retained ClinMAVE per-gene CSV exports to map curated functional assay results into `OS2`, `SBS2`, or a neutral score of `0`.

### Current Upstream Inputs

The current implementation depends on:

- `normalizedVariant.geneSymbol`
- `normalizedVariant.transcript_hgvs.mane_select_b38`
- optional `tumorType`
- local ClinMAVE files under `data/clinmave/variants.<GENE>.csv`
- local gene-role metadata in `data/_Dict_Gene.csv`

ClinMAVE `Identifier` values are normalized from forms like:

- `NM_000051.4(ATM):c.283C>T (p.Gln95Ter)`

to transcript HGVS strings like:

- `NM_000051.4:c.283C>T`

The pipeline then performs exact equality matching against `mane_select_b38`.

### Current Match Policy

- primary and only v1 match key: `normalizedVariant.transcript_hgvs.mane_select_b38`
- match requires exact string equality after ClinMAVE identifier normalization
- no fallback to alternate transcript, protein, or genomic matching in v1
- if the queried gene is not present in the retained ClinMAVE panel, return `not_available`
- if the gene is present but the variant is not found, return `not_available`

### Current Functional Class Mapping

Observed ClinMAVE functional classes in the retained dataset are:

- `Functionally normal`
- `Gain-of-function`
- `Loss-of-function`

Current interpretation:

- `Functionally normal` -> normal
- `Gain-of-function` -> GOF
- `Loss-of-function` -> LOF

### Current Gene Role Policy

Most retained genes resolve directly to oncogene or tumor suppressor gene using `data/_Dict_Gene.csv`.

`GATA3` is currently treated as a dual-role gene and requires `tumorType` to resolve role:

- oncogene contexts: `Peripheral T-Cell Lymphoma`, `T-Cell Acute Lymphoblastic Leukemia`, `Hodgkin Lymphoma`, `Neuroblastoma`, `T-Cell Lymphoblastic Lymphoma`
- tumor suppressor gene contexts: `Breast Cancer`, `Urothelial Carcinoma`, `Bladder Carcinoma`, `Renal Cell Carcinoma`, `Parathyroid Carcinoma`
- any other tumor type: no functional rule is applied and the result is `not_available`

### Current Rule Mapping

- oncogene + `Gain-of-function` -> `OS2`, score `4`
- tumor suppressor gene + `Loss-of-function` -> `OS2`, score `4`
- oncogene + `Functionally normal` -> `SBS2`, score `-4`
- tumor suppressor gene + `Functionally normal` -> `SBS2`, score `-4`
- oncogene + `Loss-of-function` -> score `0`, `applied`
- tumor suppressor gene + `Gain-of-function` -> score `0`, `applied`
- conflicting ClinMAVE classifications for the same exact matched transcript HGVS -> score `0`, `applied`

### Current Availability Semantics

- `not_available` when annotation failed
- `not_available` when no gene symbol was resolved
- `not_available` when no MANE transcript HGVS representation was resolved
- `not_available` when the gene is not in the retained ClinMAVE panel
- `not_available` when the gene is supported but the variant is not found
- `not_available` when tumor-type context is required but missing or unresolved
- `applied` when a ClinMAVE row is found and evaluated, including score `0` outcomes
- `applied` when exact-match ClinMAVE rows are found but contain conflicting classifications, in which case the statement explains that no functional rule is applied

### Current Data Access Strategy

- ClinMAVE gene CSVs are loaded lazily, one gene at a time, on first use
- parsed rows are cached in memory for the life of the process
- the pipeline does not preload all retained ClinMAVE files at startup

## OM1 Pipeline

### Purpose

The OM1 pipeline evaluates whether a localized protein-altering variant falls within a curated critical and well-established functional domain.

### Current Data Source

- local CSV: `data/om1_clingen_domains_seed.csv`
- source provenance: ClinGen CSPEC, materialized into a local MANE-anchored working table

### Current Data Access Strategy

- the service loads released OM1-ready rows from the bundled CSV on first use
- only rows with `rowStatus=ready` participate in runtime evaluation
- genes represented only by non-ready inventory rows are treated as unavailable, not neutral
- unreadable or missing local CSV data is treated as pipeline unavailability rather than as score `0`

### Current Match Policy

- matching is gene- and MANE-transcript-specific
- the pipeline requires a MANE Select transcript consequence from annotation
- eligible event types are localized protein-altering consequences currently represented as `substitution`, `deletion`, `insertion`, `duplication`, or `delins`
- substitutions use the resolved protein residue position
- localized in-frame events use the resolved protein residue span, preferring parsed `p.HGVS` event bounds and falling back to `proteinStart` and `proteinEnd`
- truncating and other non-localized event types are not eligible for the current OM1 rule and therefore evaluate to score `0` when annotation is otherwise sufficient
- excluded residues are supported by the row schema and suppress OM1 when the localized event overlaps one of those explicitly excluded positions

### Current Rule Mapping

- `OM1`, score `2`: localized MANE protein event overlaps a curated critical-domain interval in the local ClinGen table
- else score `0`

### Current Availability Semantics

- `not_available` when annotation failed
- `not_available` when the normalized gene symbol is unavailable
- `not_available` when the local OM1 domain table is missing or unreadable
- `not_available` when the gene has no released `rowStatus=ready` ClinGen domain rows in the local table
- `not_available` when no MANE Select transcript consequence is available
- `not_available` when an otherwise eligible localized protein event lacks a resolvable residue position or span
- `applied` score `0` when annotation and local OM1 data were sufficient but the current rule did not fire

## OP2 Pipeline

`OP2` is now implemented as a small curated tumor-type-aware evidence block.

### Current Rule Basis

The manuscript describes `OP2` as:

- somatic variant in a gene in a malignancy with a single genetic etiology

The current implementation intentionally treats this as a small curated exception table rather than as a broad inferred rule.

### Current Data Shape

The rules table currently lives at `data/op2_rules.csv` with columns:

- `tumorType`
- `geneSymbol`
- `maneProteinHgvs`

Current seed rows are intentionally small and representative.

### Current Match Semantics

- `tumorType` is required
- `geneSymbol` is required
- `maneProteinHgvs` is optional
- if `maneProteinHgvs` is populated, the row should match only when the MANE Select transcript consequence has the same normalized one-letter protein HGVS value
- if `maneProteinHgvs` is blank, the row acts as a broad disease-plus-gene rule

Current seeded examples include:

- variant-restricted OP2 rows for canonical disease-defining protein events such as `FOXL2 p.C134W` and `BRAF p.V600E`
- a broader disease-plus-gene OP2 row for `RB1` in retinoblastoma

### Current Availability Semantics

- `not_available` when annotation failed
- `not_available` when `tumorType` is missing because tumor context is the primary entry criterion
- `not_available` when a variant-restricted OP2 row is relevant but no MANE Select protein consequence is available for matching
- `not_available` when the local OP2 rules table is missing or unreadable
- `applied` score `1` with `OP2` when a curated row matches
- `applied` score `0` when `tumorType` is present but no curated row matches

### Current Matching Notes

- OP2 uses exact equality for `tumorType` and `geneSymbol`
- variant-restricted rows use exact equality against the MANE Select transcript consequence `proteinHgvs` value in normalized one-letter form, such as `p.V600E`
- the current broad `RB1` retinoblastoma row does not require an additional mechanism-compatibility gate

## Current Overall Score Aggregation

Current overall scoring is a three-step process:

1. Evaluate each evidence pipeline independently.
2. Apply deterministic interaction rules that suppress overlapping evidence.
3. Sum only `applied` evidence scores and map the adjusted total to a final classification.

The current evidence interaction rules are:

- suppress `OS3` if `OS1` is applicable
- suppress `OM1` if `OS1` or `OS3` is applicable
- suppress `OM4` if `OS1`, `OS3`, or `OM1` is applicable
- suppress `OM3` if `OM1` or `OM4` is applicable
- suppress `OM2` if `OVS1` is applicable

Current adjusted scoring considers these evidence lanes:

- population
- computational
- hotspots
- predictive
- om1
- op2
- functional

Suppressed and `not_available` lanes do not contribute to the final numeric score.

If all evidence lanes are `not_available`, the service does not emit a synthetic zero-score classification. Instead, the overall prediction is marked unavailable with `overallScore = null`, `overallClassification = null`, a top-level `dataAbsentReason`, and a top-level `predictionStatement`.

Current score-to-classification mapping is:

- `<= -7` -> `Benign`
- `-6..-1` -> `Likely Benign`
- `0..5` -> `VUS`
- `6..9` -> `Likely Oncogenic`
- `>= 10` -> `Oncogenic`

## Current Deferred Caveats

The manuscript's explicit exclusion rules above are implemented. Other comments and caveats from Tables 2 and 3 are currently documented as future refinement areas rather than automated logic.

More broadly, the current v1 implementation intentionally does not attempt to encode the full set of disease-specific, gene-specific, and expert-panel-specific caveats that appear across published interpretation frameworks. That includes population-threshold overrides, assay-strength downgrades, transcript and splicing caveats, domain-specific exceptions, and other special-case logic that often depends on narrow curation context. The current implementation therefore follows a smaller automated core rule set and will sometimes disagree with manually curated pipelines even when the high-level rule names appear aligned.

The highest-priority deferred caveats are:

- `OVS1` nuance for extreme 3' end pLOF variants, splice-driven in-frame rescue, alternative isoforms, and multi-transcript interpretation
- splicing-aware caution for protein-level criteria such as `OS1` and `OM4`, where the apparent amino-acid change may not reflect the true primary effect
- hotspot caution for truncating-variant-driven hotspots
- functional evidence downgrading when underlying studies are partial, conflicting, or otherwise insufficient to fully satisfy `OS2` or `SBS2`
- population-threshold refinement for hereditary cancer predisposition genes where gene-specific germline guidance should influence frequency cutoffs
- broader expert-panel and disease-specific exception handling for rules whose practical use depends on curated context beyond the generic v1 logic

## Current Prediction Rendering

The prediction endpoints currently return one FHIR Observation-style object per variant.

Current rendering intent:

- `Observation.code`: temporary code for oncogenicity prediction
- `Observation.issued`: timestamp when the service generated the prediction
- `Observation.extension`: custom extension carrying the originally submitted variant string in `valueString`
- `Observation.valueInteger`: overall numeric score when available
- `Observation.interpretation`: final overall classification when available
- `Observation.dataAbsentReason`: top-level unavailable-prediction signal when no evidence lanes are evaluable
- zero or more `Observation.component` entries for score-contributing evidence
- `component.code`: temporary code identifying the pipeline
- `component.valueInteger`: pipeline score for an included applied evidence lane
- `component.interpretation.coding.code`: pipeline evidence code such as `OP4`, `SBS1`, or `SBVS1`
- `component.interpretation.text`: short clinician-facing evidence statement

The FHIR projection is intentionally compact:

- include only evidence components where `status == "applied"` and `score != 0`
- omit `suppressed` lanes, `not_available` lanes, and `applied` neutral lanes with score `0`
- when no evidence lanes are evaluable, omit `valueInteger` and use top-level `dataAbsentReason` instead of forcing a numeric score/classification
- rely on `GET /summarizeEvidence` for the full audit surface, including suppression and availability details

Current serialization behavior:

- prediction endpoints omit `null` fields from the JSON response
- defaulted fields such as `resourceType="Observation"` and `status="final"` are still emitted

The final client-facing result should not embed `AnnotatedVariant` or `NormalizedVariant`. Any normalization or annotation provenance needed by the client should be surfaced through the evidence output itself, primarily through fields such as:

- `source`
- `matchedData`
- `evidenceStatement`

## Known v1 Gaps

- some evidence policies remain intentionally narrow, especially exact-match functional lookups and the small curated OP2 rule table
- many disease-specific, gene-specific, and expert-panel-specific caveats remain documented limitations rather than automated logic
- generic thresholds and rule mappings are still used in places where mature frameworks apply narrower population, transcript, domain, or assay-specific exceptions

In the current REST response shape, the FATHMM-family fields exposed for this implementation are the `FATHMM-XF` dbNSFP keys with hyphenated names such as `fathmm-xf_coding_pred`.
Although those field names are the ones we read from the response, the current Ensembl REST service returned `invalid_field` when they were requested explicitly, so the implementation uses `dbNSFP=ALL` and then extracts the needed keys from the response.

### Current Intended Rule Basis

- `OP1` is used as a positive supporting rule when `CADD PHRED >= 15`
- `SBP1` is used as a benign supporting rule when both of the following are true:
- `mostSevereConsequence == "missense_variant"`
- `CADD PHRED < 15`
- `FATHMM-XF` prediction is benign or neutral

Low `CADD` alone does not trigger `SBP1`. If `CADD` is low but the `FATHMM-XF` call is missing or not benign/neutral, the pipeline returns a score of `0` with no evidence code. Non-missense variants can still receive `OP1` if `CADD` is high enough, but they are not eligible for the current `SBP1` benign rule.

### Initial Rule Mapping

- If annotation failed or computational inputs are unavailable, return `not_available` with score `0`
- If `CADD PHRED >= 15`, return `OP1` with score `1`
- Else if `mostSevereConsequence == "missense_variant"` and `FATHMM-XF` is benign or neutral, return `SBP1` with score `-1`
- Else if `mostSevereConsequence` is not `missense_variant`, return score `0` with no evidence code
- Else return score `0` with no evidence code

### Draft Decision Logic

```text
if annotation failed:
    not_available
elif cadd_phred is missing:
    not_available
elif cadd_phred >= 15:
    score = 1
    evidenceCode = "OP1"
elif most_severe_consequence == "missense_variant" and fathmm_xf_prediction in {"N", "neutral", "benign", "tolerated"}:
    score = -1
    evidenceCode = "SBP1"
elif most_severe_consequence != "missense_variant":
  score = 0
  evidenceCode = null
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
  "evidenceStatement": "CADD supports oncogenicity for this variant (PHRED 25.3; most severe consequence missense_variant).",
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
  "evidenceStatement": "Computational missense benign rules were not applicable because the most severe consequence was synonymous_variant.",
  "status": "applied",
  "source": "vep"
}
```

## Hotspots Pipeline

Planned source: the local workbook at `data/hotspots_v2.xlsx`, copied from Cancer Hotspots: https://www.cancerhotspots.org/#/home

Current implementation:

- the workbook is loaded into an in-memory cache on first use rather than being reopened for each variant
- transcript consequences are evaluated in order, and the first defensible hotspot match is used
- matching is gene-centric because the workbook is not transcript-indexed
- SNVs require exact match on gene, amino-acid position, reference amino acid, and alternate amino acid
- indels are supported only for direct matches that can be justified from the workbook's native event representation, with no additional indel normalization layer
- thresholds follow the legacy implementation:
- `OS3`: `mutation_count >= 50` and exact protein-event count `>= 10`
- `OM3`: exact protein-event count `>= 10`
- `OP3`: exact protein-event count `1..9`

## Predictive Pipeline

The predictive pipeline is currently implemented and combines deterministic consequence-based logic with ClinVar somatic oncogenicity lookups.

Current rule order:

- `OVS1`
- `OS1`
- `OM2`
- `SBP2`
- `OM4`
- else score `0` with `applied`

Current source behavior:

- consequence-driven rules use normalized and annotated variant state
- role-sensitive consequence rules use optional `tumorType` to resolve dual-role genes such as `GATA3`; without a resolved role, `OVS1` is not applied
- ClinVar lookup uses E-utilities search plus summary retrieval
- ClinVar `[varnam]` matching is treated as non-exact, so local alias confirmation is required before applying `OS1` or `OM4`

## Functional Pipeline

### Purpose

The functional pipeline uses locally retained ClinMAVE per-gene CSV exports to map curated functional assay results into `OS2`, `SBS2`, or a neutral score of `0`.

### Current Upstream Inputs

The current implementation depends on:

- `normalizedVariant.geneSymbol`
- `normalizedVariant.transcript_hgvs.mane_select_b38`
- optional `tumorType` request input on evidence and prediction endpoints
- local ClinMAVE files under `data/clinmave/variants.<GENE>.csv`
- local gene-role metadata in `data/_Dict_Gene.csv`

ClinMAVE `Identifier` values are normalized from forms like:

- `NM_000051.4(ATM):c.283C>T (p.Gln95Ter)`

to transcript HGVS strings like:

- `NM_000051.4:c.283C>T`

The pipeline then performs exact equality matching against `mane_select_b38`.

### Current Match Policy

- Primary and only v1 match key: `normalizedVariant.transcript_hgvs.mane_select_b38`
- Match requires exact string equality after ClinMAVE `Identifier` normalization
- No fallback to alternate transcript, protein, or genomic matching in v1
- If the queried gene is not present in the retained ClinMAVE panel, return `not_available`
- If the gene is present but the variant is not found, return `not_available`

### Current Functional Classification Mapping

Observed ClinMAVE functional classes in the retained dataset are:

- `Functionally normal`
- `Gain-of-function`
- `Loss-of-function`

They are interpreted as:

- `Functionally normal` -> normal
- `Gain-of-function` -> GOF
- `Loss-of-function` -> LOF

### Current Gene Role Policy

Most retained genes resolve directly to `oncogene` or `tsg` using `data/_Dict_Gene.csv`.

`GATA3` is currently treated as a dual-role gene and requires `tumorType` to resolve role:

- Oncogene contexts: `Peripheral T-Cell Lymphoma`, `T-Cell Acute Lymphoblastic Leukemia`, `Hodgkin Lymphoma`, `Neuroblastoma`, `T-Cell Lymphoblastic Lymphoma`
- Tumor suppressor contexts: `Breast Cancer`, `Urothelial Carcinoma`, `Bladder Carcinoma`, `Renal Cell Carcinoma`, `Parathyroid Carcinoma`
- Any other tumor type: no functional rule is applied and the result is `not_available`

### Current Rule Mapping

- oncogene + `Gain-of-function` -> `OS2`, score `4`
- tumor suppressor gene + `Loss-of-function` -> `OS2`, score `4`
- oncogene + `Functionally normal` -> `SBS2`, score `-4`
- tumor suppressor gene + `Functionally normal` -> `SBS2`, score `-4`
- oncogene + `Loss-of-function` -> score `0`, `applied`
- tumor suppressor gene + `Gain-of-function` -> score `0`, `applied`
- conflicting ClinMAVE classifications for the same exact matched transcript HGVS -> score `0`, `applied`

### Current Availability Semantics

- `applied` when a ClinMAVE row is found and evaluated, including neutral score `0` outcomes
- `applied` when exact-match ClinMAVE rows are found but contain conflicting functional classifications, in which case the statement explains that no functional rule is applied
- `not_available` when the gene is not in the retained ClinMAVE panel
- `not_available` when the gene is supported but the variant is not found
- `not_available` when tumor-type context is required but missing or unresolved

### Current Data Access Strategy

- ClinMAVE gene CSVs are loaded lazily, one gene at a time, on first use
- Parsed rows are cached in memory for the life of the process
- The pipeline does not preload all retained ClinMAVE files at startup

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