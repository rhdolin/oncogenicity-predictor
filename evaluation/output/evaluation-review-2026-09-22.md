# Evaluation Review By Reference Set

Date: 2026-09-22

This review covers the three suffix groups currently present in `evaluation/output/`:

- `ClinVar`
- `Horak`
- `OncoV1Manual`

The goal here is not just to restate the metrics, but to interpret what they imply about algorithm behavior, likely caveats, and revision priorities.

## Executive Summary

Three broad patterns are consistent across the reference sets.

Before comparing the groups directly, one caveat needs to be stated up front: the ClinVar evaluation is not fully independent because ClinVar-derived information is also used by the algorithm. That creates a circularity risk and likely makes the ClinVar results look more favorable than they would under a completely external reference set. ClinVar is still useful for stress-testing behavior, but it should not be treated as the strongest evidence of generalization.

First, the algorithm is generally conservative. Across all three groups, the dominant directional error is undercalling oncogenic variants into `VUS`, not overcalling benign or `VUS` variants into oncogenic categories.

Second, extreme benign-versus-oncogenic reversals are uncommon. That is good news. The main problem is not gross polarity inversion across the board; it is failure to assign enough positive weight to many reference-set oncogenic variants.

Third, where score-level evaluation is available, the top-line category agreement is materially better than the score concordance. That means the pipeline often lands in the right neighborhood but does not reproduce the reference-set scoring rationale cleanly. Exact criteria-set agreement is especially low.

The most important cross-group algorithm revision themes are:

- reduce the tendency to collapse positive variants into `VUS`
- audit benign evidence rules that appear to overpower known positive evidence in some genes
- improve functional and predictive swimlane calibration
- decide whether the intended scope includes broad loss-of-function interpretation outside the current core driver-gene logic
- fix annotation failures that currently block a nontrivial portion of the ClinVar run

## Top-Line Comparison

| Reference set | Rows | Unavailable | Observed agreement | Weighted agreement linear | Weighted kappa linear | Adjacent disagreements | Extreme disagreements | Exact score agreement | Score > 2 off |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ClinVar | 316 | 20 | 0.882 | 0.936 | 0.634 | 32 | 3 | n/a | n/a |
| Horak | 94 | 0 | 0.713 | 0.846 | 0.566 | 25 | 2 | 0.149 | 0.447 |
| OncoV1Manual | 135 | 0 | 0.719 | 0.856 | 0.315 | 37 | 1 | 0.163 | 0.348 |

Important context for comparing these rows:

- The ClinVar group is category-only, so it supports final classification review but not score-level concordance.
- The ClinVar group is also partially self-referential because ClinVar-derived information contributes to the predictive evidence framework. Its concordance is therefore likely somewhat optimistic relative to a fully external reference set.
- The Horak and OncoV1Manual groups support full score-based review, but they were derived differently and should not be treated as interchangeable benchmarks.
- The OncoV1Manual group has an additional scoring-framework mismatch: its reference-set scoring can include multiple criteria that map to the same conceptual swimlane, whereas the current algorithm emits at most one criterion per swimlane. That makes direct score concordance harsher for OncoV1Manual than for a reference set built around the same one-criterion-per-swimlane architecture.
- The much lower kappa in OncoV1Manual does not simply mean the model is worse there. It also reflects the different category distribution in that reference set, especially the heavy `VUS` representation.

## Cross-Group Findings

### 1. Undercalling is the dominant systematic error

The main recurring issue is not that the algorithm invents oncogenicity too freely. It is that it often fails to accumulate enough positive evidence to preserve an oncogenic reference-set label.

Evidence for that pattern:

- ClinVar: among evaluable oncogenic-family reference rows, `226` remained in the oncogenic family, but `22` dropped to `VUS` and `3` dropped to benign.
- Horak: among oncogenic-family reference rows, only `21` stayed oncogenic-family, while `22` dropped to `VUS` and `2` to benign.
- OncoV1Manual: among oncogenic reference rows, only `16` stayed oncogenic, `28` dropped to `VUS`, and `1` dropped to benign.

This is the single clearest algorithmic trend across the runs.

### 2. The final category often looks better than the scoring rationale underneath it

For the score-evaluable sets, exact criteria-set agreement is low:

- Horak: `5.3%`
- OncoV1Manual: `10.4%`

That means the algorithm often arrives near the reference-set conclusion without reproducing the same evidence pattern. In practice, that creates two risks:

- the model may look more aligned than it really is if only final classification is examined
- future changes could destabilize apparently good category-level performance because the underlying rationale is already not tightly anchored

### 3. Population and `op2` are stable; functional and predictive are the main disagreement sources

The score-based swimlanes show a consistent structure.

Population is already very strong.

- Horak: exact `87.2%`, acceptable `95.7%`, large deviation `4.3%`
- OncoV1Manual: exact `99.3%`, acceptable `100%`, large deviation `0%`

`op2` is perfect in both score-evaluable sets.

Functional is the weakest swimlane.

- Horak: exact `57.4%`, acceptable `57.4%`, large deviation `42.6%`
- OncoV1Manual: exact `79.3%`, acceptable `79.3%`, large deviation `20.7%`

The directionality matters here. Functional disagreement is mostly downward rather than upward. In Horak, the functional swimlane has a mean predicted-minus-reference delta of about `-1.36`, with `35` downward versus `5` upward disagreements. In OncoV1Manual, the same pattern persists more mildly, with mean delta about `-0.33` and `19` downward versus `9` upward disagreements. So functional is not just noisy; it is usually underscoring relative to the reference sets.

Predictive is the next most important disagreement source.

- Horak: large deviation `16.0%`
- OncoV1Manual: large deviation `15.6%`

Predictive disagreement is less uniform across reference sets. In Horak it trends upward slightly overall, with mean delta about `+0.44` and more upward than downward disagreements, suggesting that predictive is sometimes compensating for missing evidence elsewhere. In OncoV1Manual it trends downward, with mean delta about `-0.86` and `21` downward versus `7` upward disagreements, which is more consistent with the broader undercalling pattern in that reference set.

Hotspots is moderate in Horak and relatively strong in OncoV1Manual.

### 4. Benign-evidence rules need a targeted audit

Several of the most concerning disagreements involve benign evidence appearing to overpower known positive context.

Examples include:

- TP53 variants in ClinVar and Horak that are downgraded toward `Likely Benign`
- a BRCA2 variant in OncoV1Manual downgraded from oncogenic to `Likely Benign`
- an RB1 variant in Horak downgraded from `Likely Oncogenic` to `Benign`

In multiple outliers, the emitted criteria include `SBS2` or `SBVS1`. Those rules deserve direct review to confirm that their gating logic is not too permissive in somatic-driver contexts.

More concretely, the two rules do not look equally problematic.

#### `SBS2`: a recurring source of discordance, but not necessarily a rule defect

`SBS2` is the more important benign-side discordance source, but that does not automatically mean the rule is wrong.

Observed pattern:

- ClinVar: `SBS2` appears in `21` predictions, `19` of which are category mismatches, including `7` large disagreements
- Horak: `SBS2` appears in `9` predictions, `5` of which are mismatches, including `2` large disagreements
- OncoV1Manual: `SBS2` appears only once, and that one case is a large disagreement

Representative severe examples:

- ClinVar TP53 variants downgraded from oncogenic-family labels to `VUS` or `Likely Benign` with `SBS2` present
- ClinVar PTEN variants downgraded from `Oncogenic` to `VUS` with `OS1` and `SBS2` both present
- Horak `NM_000546.5:c.1040C>A` in TP53: reference `Likely Oncogenic`, predicted `Likely Benign` with `OP4 + OP1 + SBS2`
- Horak `NM_004333.6:c.1780G>A` in BRAF: reference `Oncogenic`, predicted `VUS` with `OP4 + OP1 + OS3 + SBS2`
- OncoV1Manual `NC_000013.11:g.32362681A>G` in BRCA2: reference `Oncogenic`, predicted `Likely Benign` with `OP4 + OP1 + SBS2`

What the current code is doing:

- `SBS2` is applied whenever ClinMAVE finds an exact MANE transcript match and the functional classification is `Functionally normal` in a resolved oncogene or tumor suppressor gene context.
- It always contributes `-4` points.
- It is not downweighted for coexistence with strong positive evidence.
- It is not suppressed downstream. The current score calculator suppresses some positive redundancies (`OS3`, `OM1`, `OM4`, `OM2` in selected contexts), but it never suppresses `SBS2`.

This does not, by itself, prove a bug. If one accepts the scoring premise that a high-quality functional result showing normal activity should contribute `-4`, then at least some of these disagreements may reflect legitimate evidence conflicts between ClinMave-backed functional evidence and the reference-set interpretation rather than an implementation error.

The concrete review target is therefore narrower than "fix `SBS2`". The review should ask:

1. Was the ClinMave match correct for the exact variant and transcript?
2. Is the underlying functional result strong enough and relevant enough to warrant a full `SBS2` assignment?
3. Is the reference-set classification intentionally outweighing normal-function evidence with other somatic signals?
4. Are these cases examples of real biological tension, where normal activity in one assay does not fully negate somatic oncogenic relevance?

If those questions are answered satisfactorily, then the disagreement should be documented as a reference-set conflict rather than treated as an algorithm bug. Only if the match, assay applicability, or use of the rule is incorrect does `SBS2` become a concrete algorithm revision target.

#### `SBVS1`: probably not a broad rule failure, but still worth a conflict audit

`SBVS1` appears to be a smaller problem.

Observed pattern:

- ClinVar: `SBVS1` appears in `3` predictions, with `1` mismatch and no large disagreements
- Horak: `SBVS1` appears in `5` predictions, with `1` mismatch and `1` large disagreement
- OncoV1Manual: `SBVS1` does not appear in the evaluated predictions

So unlike `SBS2`, `SBVS1` is not showing up repeatedly as a systematic source of discordance.

The more important caveat is about how the rule is triggered. The current implementation uses the maximum subpopulation allele frequency when available, and only falls back to the overall allele frequency if no subpopulation value is present. That means `SBVS1` can fire because a variant is common in one ancestry-specific gnomAD bucket even when the overall allele frequency is much lower. This is not necessarily wrong, but it is an interpretation choice that should be stated explicitly because it may surprise readers who assume the rule is based on overall gnomAD frequency alone.

What the current code is doing:

- it takes a single effective allele frequency, preferring max subpopulation AF over overall AF
- if that value is `> 5%`, it assigns `SBVS1` with score `-8`
- it does not currently consider gene role, somatic recurrence, hotspot status, or coexistence with strong positive evidence
- it is not suppressed downstream when strong positive evidence is present

One concrete example of the aggregation-policy caveat is the Horak PTEN case, where `SBVS1` fires from a max subpopulation allele frequency of about `5.33%` even though the overall allele frequency is only about `0.20%`. Again, that may be the intended policy, but it should be transparent in the write-up.

That said, the observed failure pattern suggests that `SBVS1` is usually not wrong by itself. The bigger issue is that when it fires, the rest of the pipeline may fail to recover enough positive evidence to offset it. In other words, `SBVS1` looks more like a conflict-resolution problem than the main driver of error.

The concrete review items for `SBVS1` should therefore be:

1. Confirm that the high-frequency calls behind `SBVS1` are biologically plausible and not annotation artifacts.
2. Decide whether the intended `SBVS1` policy is truly "common in any sufficiently represented subpopulation" rather than "common overall," and document that choice explicitly.
3. Add an explicit conflict check for cases where `SBVS1` coexists with very strong positive somatic evidence, so those rows are surfaced for manual review or special handling rather than blindly summed.
4. Decide whether `SBVS1` should remain fully counted when recurrent somatic-driver evidence is present, or whether such combinations should be classified as internally inconsistent and handled separately.

The most likely conclusion here is that `SBVS1` should remain available, but the algorithm needs a clearer policy for what to do when a very strong benign population signal and a very strong somatic signal coexist.

## Reference Set Review: ClinVar

### Independence caveat

The ClinVar group is useful, but it is not a fully independent evaluation. ClinVar-derived information contributes to the algorithm's predictive evidence, so some variants in this reference set can receive favorable support from a source that is also represented in the reference labels.

That creates a circularity risk. In practical terms, the ClinVar metrics likely overestimate external performance, especially for variants that benefit from ClinVar-informed predictive evidence such as prior oncogenic assertions or same-amino-acid-change matching. The ClinVar results should therefore be interpreted as partly self-referential and should not be weighted as heavily as an independent external reference set when summarizing overall performance.

### Overall reading

The ClinVar run is the strongest of the three on top-line category concordance. Observed agreement is `88.2%`, weighted linear agreement is `93.6%`, and weighted quadratic agreement is `96.3%`. That means the model is usually in the correct or nearly correct ordinal neighborhood relative to this reference set.

However, the kappa values are more modest at `0.63` to `0.65`, which means some of the apparent agreement is helped by the class distribution. More importantly, because the reference set is not independent of one of the evidence sources, this is good performance but not as definitive as the raw agreement alone may suggest.

### What the disagreements mean

This run is best described as family-level concordant but boundary-unstable.

The most common disagreement pattern is swapping `Oncogenic` and `Likely Oncogenic`, not flipping all the way across the scale. That means the algorithm often recognizes oncogenicity, but is inconsistent about whether the evidence is strong enough to cross from likely to definite oncogenic.

There is also a secondary undercalling pattern where oncogenic-family reference rows fall to `VUS`. That is less frequent than in Horak or OncoV1Manual, but it is still present and clinically more important than the likely-versus-definite boundary swaps.

### Operational caveat

The main ClinVar-specific caveat is not scoring. It is availability.

There are `20` unavailable predictions out of `316` rows, or `6.3%` of the batch. All of them appear to be annotation failures with `dataAbsentReason = error`.

The failures are concentrated in a small set of genes:

- TP53: `7`
- PTEN: `6`
- PIK3CA: `3`
- BRAF: `2`

That concentration strongly suggests an implementation issue rather than random noise.

### Specific concerns

The most concerning category reversals in this batch are TP53 downgrades from oncogenic-family labels to benign-family labels. Those are rare, but they should be treated as priority review cases because they are exactly the sort of errors that can undermine confidence even when aggregate performance is otherwise good.

Several additional exact large disagreements also involve TP53, PTEN, and FLT3. That concentration suggests that gene-specific logic is more likely involved than broad random miscalibration.

### Review conclusion for ClinVar

This reference set suggests that the pipeline is directionally sound at the final classification level, but still too conservative at some oncogenic boundaries and not operationally robust enough yet for certain variant classes or loci. If the ClinVar run were the only evidence, the biggest concerns would be:

- fixing annotation failures
- reducing the small but important set of severe TP53-like downgrades
- tightening the distinction between `Likely Oncogenic` and `Oncogenic`

## Reference Set Review: Horak

### Overall reading

The Horak run is materially weaker than ClinVar and exposes the underlying scoring problems more clearly.

Observed agreement is `71.3%`, weighted linear agreement is `84.6%`, and weighted kappa linear is `0.57`. That profile says the model is often near the reference-set category, but not at a level that should be described as strong ordinal alignment by the stated `>= 0.90` weighted-linear goal.

### Main behavior: strong positive undercalling

The clearest message from Horak is that the algorithm undercalls oncogenic variants.

Among oncogenic-family reference rows:

- `21` remain oncogenic-family
- `22` are downgraded to `VUS`
- `2` are downgraded to benign-family labels

By contrast, `VUS` is comparatively stable and benign-family reference rows are mostly handled well. So the algorithm is not broadly confused; it is specifically reluctant to preserve higher-confidence positive calls.

### Score-level interpretation

The score metrics reinforce that conclusion.

- Exact score agreement is only `14.9%`.
- Agreement within 2 points is `55.3%`.
- Large deviations greater than 2 points occur in `44.7%` of rows.

That is too high. It means nearly half of the scored cases are not merely slightly off; they are meaningfully misaligned with the reference set.

The exact criteria-set agreement rate is only `5.3%`, which is another sign that the model is not reproducing the reference-set reasoning structure.

### Swimlane interpretation

Population is strong and not a major concern.

Computational and `om1` have no large deviations, but they are somewhat fuzzy rather than consistently exact. They are serviceable, not the main bottleneck.

The biggest problems are:

- functional: `42.6%` large deviations
- predictive: `16.0%` large deviations
- hotspots: `12.8%` large deviations

Functional is the clearest weakness in this reference set, and it is predominantly a down-scoring problem rather than bidirectional noise. The functional swimlane has a mean predicted-minus-reference delta of about `-1.36`, with `35` downward disagreements and only `5` upward disagreements. In other words, when functional misses in Horak, it usually withholds positive weight that the reference set expected.

Hotspots also trends downward in Horak, not upward. Its mean delta is about `-0.62`, with `19` downward disagreements and no upward disagreements. That suggests the hotspot logic is more often failing to credit recurrent context than inventing it.

Predictive behaves differently. Its large-deviation rate is still meaningful, but its net direction is slightly upward in Horak, with mean delta about `+0.44`. That means predictive is not the main source of undercalling here; if anything, it partly offsets shortfalls elsewhere.

### Concrete disagreement pattern

The largest negative score deltas in Horak are driven by variants where the reference set assigns strong positive evidence and the model emits only weak or moderate support.

Examples include:

- TP53 `NM_000546.5:c.1040C>A`: reference `Likely Oncogenic` score `7`, predicted `Likely Benign` score `-2`
- BRAF `NM_004333.6:c.1780G>A`: reference `Oncogenic` score `10`, predicted `VUS` score `2`
- TERT promoter variants: reference `Likely Oncogenic` score `9`, predicted `VUS` score `1`
- FLT3 variants: reference `Oncogenic` score `10`, predicted `VUS` score `3`
- EZH2 variants: reference `Likely Oncogenic` or `Oncogenic`, predicted `VUS`

Those are not random misses. They suggest missing or suppressed strong positive rules in well-known recurrent contexts.

### Review conclusion for Horak

The Horak run suggests a real algorithmic calibration problem, not just cosmetic disagreement. The strongest evidence points to insufficient positive weighting in functional and predictive interpretation, plus occasional over-application of benign evidence. This run should be treated as a primary driver for algorithm revision.

## Reference Set Review: OncoV1Manual

### Overall reading

The OncoV1Manual run has similar raw category agreement to Horak, but a much lower kappa.

- Observed agreement: `71.9%`
- Weighted linear agreement: `85.6%`
- Weighted kappa linear: `0.315`

The low kappa means that once the category distribution is taken into account, the apparent agreement is not especially strong. This likely reflects both algorithm behavior and reference-set composition.

There is also a structural comparability caveat in this reference set: the reference-set scoring is not constrained by the same one-criterion-per-swimlane rule that the current algorithm uses. The current implementation emits at most one criterion from each swimlane, while OncoV1Manual can include multiple criteria that would collapse into the same swimlane in this repository.

For example, the XRCC2 row `NC_000007.14:g.152660728del` is labeled `Oncogenic` with score `18` and criteria `[`OVS1`, `OS1`, `OS2`, `OP1`, `OP4`]`. In the current algorithm, `OVS1` and `OS1` both live in the predictive swimlane and cannot both contribute simultaneously. So some of the score gap in this reference set is not just evidence disagreement; it is a framework mismatch between the reference-set scoring model and the implementation being evaluated.

### Main behavior: VUS inflation

This run is dominated by a `VUS` sink.

Reference-set distribution:

- `86` VUS
- `45` oncogenic
- `4` benign

Predicted distribution:

- `113` VUS
- `21` oncogenic-family
- `1` benign-family

Among oncogenic reference rows:

- `16` remain oncogenic-family
- `28` drop to `VUS`
- `1` drops to benign-family

This is the strongest undercalling signal of the three groups.

### Score-level interpretation

The score-level picture is slightly better than Horak, but still not good enough to call robust.

- Exact score agreement: `16.3%`
- Within 2 points: `65.2%`
- Large deviations greater than 2 points: `34.8%`

Exact criteria-set agreement improves slightly to `10.4%`, but remains low.

That low exact criteria-set agreement should be interpreted with caution here, because OncoV1Manual is at least partly asking the current implementation to reproduce combinations of criteria that it is architecturally unable to emit. In other words, some of the score-level discordance is true evidence disagreement, but some of it is a consequence of comparing two different scoring grammars.

### Swimlane interpretation

Population is excellent.

Hotspots is also comparatively strong in this set, with only `1.5%` large deviations.

Computational and `om1` again look acceptable but fuzzy rather than exact.

The main pressure points are:

- functional: `20.7%` large deviations
- predictive: `15.6%` large deviations

Again, directionality is important. Functional remains net downward, though less severely than in Horak. Its mean predicted-minus-reference delta is about `-0.33`, with `19` downward and `9` upward disagreements. Predictive is more clearly downward in this reference set, with mean delta about `-0.86` and `21` downward versus `7` upward disagreements. That makes predictive, together with `om1`, a plausible driver of the strong `VUS` sink seen in this batch.

`om1` also deserves explicit mention here. Although it is not one of the highest large-deviation swimlanes by the current summary metric, it is directionally one-sided in OncoV1Manual: mean delta about `-0.86`, with `58` downward and `0` upward disagreements. That pattern is consistent with the possibility that the algorithm's `om1` usage is systematically narrower than the reference-set interpretation model.

### Concrete disagreement pattern

The most informative feature of this run is the outlier list. Several of the largest negative score deltas involve loss-of-function or apparently high-impact variants in genes such as:

- XRCC2
- BARD1
- NBN
- FANCL
- EPHA7
- ANKRD11
- BRIP1

In many of these rows, the reference set assigns strong positive evidence such as `OVS1`, `OS2`, or `OM1`, while the model emits only `OP4` or `OP4 + OP1` and leaves the case at `VUS`.

That can mean one of two things:

- the algorithm is missing intended loss-of-function positive logic
- the reference set includes biology that is broader than the current algorithm scope

That distinction matters. This is the reference set where scope clarification is most important.

### Review conclusion for OncoV1Manual

This run does not merely show generic noise. It shows a structural mismatch between the model’s current positive-evidence logic and the kinds of variants that this reference set treats as oncogenic. If this reference set is in scope, the algorithm needs a broader and stronger framework for loss-of-function and gene-role-aware positive evidence. If it is partly out of scope, that limitation should be stated explicitly.

It also shows a scoring-framework mismatch that should be stated plainly: some OncoV1Manual rows use multiple positive criteria that would occupy the same swimlane in the current implementation, so perfect score replication is not achievable without changing the model architecture.

## Suggested Algorithm Revisions

### 1. Fix annotation failures before tuning concordance claims

ClinVar currently loses `6.3%` of rows to annotation failure. That is too high to treat as a minor operational detail. These failures cluster in specific genes and should be addressed with targeted debugging, retries, and regression tests.

### 2. Reduce the `VUS` sink for oncogenic reference-set variants

Across all groups, the most common clinically important miss is `Oncogenic` or `Likely Oncogenic` collapsing to `VUS`. That points to insufficient accumulation of strong positive evidence, not rampant false-positive behavior.

### 3. Audit `SBS2` and `SBVS1` gating logic

Several severe downgrades include these benign rules. Review whether they should be:

- disallowed in certain recurrent somatic-driver contexts
- weakened when strong positive hotspot or predictive evidence is present
- made more gene-role aware

### 4. Prioritize functional swimlane revision

Functional is the weakest score-evaluable swimlane in both Horak and OncoV1Manual. It is the clearest contributor to large score deviations and should be treated as a first-line revision target.

### 5. Strengthen predictive logic for known recurrent contexts

The Horak outliers strongly suggest missing or underweighted positive logic for cases such as TERT promoter events, recurrent BRAF and FLT3 alterations, EZH2 hotspot variants, and selected TP53 contexts.

### 6. Decide the intended scope for loss-of-function interpretation

The OncoV1Manual outliers imply that the current algorithm may be narrower than the reference-set interpretation model. If broad loss-of-function oncogenic inference is intended, then `OVS1`-like and related gene-role-aware logic needs expansion. If it is not intended, the limitation should be documented clearly.

### 7. Add a targeted regression panel

A good next step would be a small hand-curated regression panel containing:

- ClinVar unavailable variants
- TP53 benign-family reversals
- Horak high-negative-delta variants such as TERT, FLT3, EZH2, and BRAF
- OncoV1Manual loss-of-function outliers such as XRCC2, BARD1, NBN, FANCL, and BRCA2

That would let you test revisions against the exact failures driving the current review.

## Final Assessment

The current implementation already shows useful signal. It is not randomly classifying variants, and it rarely produces catastrophic benign-versus-oncogenic reversals. That is the strongest positive takeaway.

The main limitation is conservative undercalling, especially when the reference set expects stronger positive evidence than the current rules are emitting. ClinVar looks closest to operational usefulness at the final category level, but that result is partly confounded by source overlap and it still has an availability problem. Horak and OncoV1Manual show that the underlying scoring logic is not yet tightly aligned with the reference-set reasoning, especially in functional and predictive evidence.

So the right near-term message is:

- category-level behavior is directionally promising
- score-level concordance is still incomplete
- functional and predictive evidence are the main calibration targets
- gene-context-specific downgrades and annotation failures should be treated as the highest-priority fixes