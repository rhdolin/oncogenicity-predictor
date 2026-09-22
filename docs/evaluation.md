# Evaluation Workflow

## Purpose

The local evaluation workflow compares the current FHIR prediction output against a curated reference set.

For v1, evaluation is intentionally based on the final surfaced FHIR criteria, overall score, and overall classification rather than on the full internal reasoning state. The full internal summary is still preserved in the output table for deep-dive review.

## Input Files

- `evaluation/variantLists/variantList.csv`: selected reference set used for routine evaluation
- `evaluation/variantLists/variantListMaster.csv`: larger source pool from which `variantList.csv` can be curated

Both files retain the legacy field names:

- `variant`
- `gene`
- `source`
- `classification`
- `points`
- `criteria`
- `comments`

Both files also include:

- `referenceDetailLevel`

`evaluation/variantLists/variantList.csv` is the canonical evaluator input and supports two row types:

- fully specified rows with integer `points` values and JSON-list `criteria` values
- category-only rows where `classification` is present, `points` and `criteria` are blank, and any legacy free-text criteria note is preserved in `comments`

`evaluation/variantLists/variantListMaster.csv` uses the same mixed-detail schema. Fully specified rows are normalized to integer `points` values plus JSON-list `criteria` values, while rows sourced from category-only references are marked with `referenceDetailLevel = category_only`, leave `points` and `criteria` blank, and preserve the original free-text criteria note in `comments`.

## Runner

Run the evaluation locally from the repository root:

```bash
source .venv/bin/activate
python evaluation/runEvaluation.py
```

The runner:

1. Loads `evaluation/variantLists/variantList.csv`.
2. Runs the local normalization, annotation, evidence, scoring, and FHIR projection pipeline for each row.
3. Stores the exact `/summarizeEvidence`-style JSON payload in the `evidenceSummary` column.
4. Flattens the final FHIR Observation into one row per input variant.
5. Overwrites the stable output files in `evaluation/output/`.

## Output Files

### `evaluation/output/oncogenicityPredictions.csv`

One row per input variant with these columns:

- `variant`
- `gene`
- `source`
- `referenceClassification`
- `referenceScore`
- `referenceCriteria`
- `predictedClassification`
- `predictedScore`
- `predictedCriteria`
- `predictionUnavailable`
- `evidenceSummary`

Column notes:

- `referenceCriteria` and `predictedCriteria` are serialized as JSON lists in a single cell.
- `predictionUnavailable` is a boolean serialized as JSON `true` or `false`.
- `evidenceSummary` stores the full internal summary JSON used to build the FHIR projection.
- For category-only reference rows, `referenceScore` is blank and `referenceCriteria` is `[]`.

### `evaluation/output/metrics-category-based.csv`

This file summarizes overall classification concordance.

Current contents include:

- input-row count
- unavailable-prediction count
- category-tested row count
- score-tested row count
- observed agreement
- weighted agreement and weighted kappa
- adjacent and extreme disagreement counts
- exact overall score agreement
- overall score agreement within 2 points
- a 3x3 confusion matrix using reference buckets `Benign`, `VUS`, and `Oncogenic`

Category-only reference rows remain part of the classification confusion matrix and agreement statistics. They are excluded only from the score-specific calculations.

Field definitions:

- `numInputRows`: total number of input rows processed by the evaluator
- `numPredictionUnavailable`: number of rows where no final prediction was returned and the row was excluded from concordance metrics
- `numCategoryTested`: number of rows included in category-based comparison
- `numScoreTested`: number of rows included in score-based comparison
- `observedAgreement`: exact category agreement rate, computed as `diagonalAgreement / numCategoryTested`
- `weightedAgreementLinear`: ordinal agreement with partial credit for adjacent disagreements, computed as `(diagonalAgreement + 0.5 * adjacentDisagreements) / numCategoryTested`
- `expectedWeightedLinear`: chance-expected ordinal agreement using the same linear weighting scheme
- `weightedKappaLinear`: chance-corrected ordinal agreement using linear weights, computed as `(weightedAgreementLinear - expectedWeightedLinear) / (1 - expectedWeightedLinear)`
- `weightedAgreementQuadratic`: ordinal agreement with stronger partial credit for near-misses, computed as `(diagonalAgreement + 0.75 * adjacentDisagreements) / numCategoryTested`
- `expectedWeightedQuadratic`: chance-expected ordinal agreement using quadratic-style weights
- `weightedKappaQuadratic`: chance-corrected ordinal agreement using quadratic-style weights, computed as `(weightedAgreementQuadratic - expectedWeightedQuadratic) / (1 - expectedWeightedQuadratic)`
- `adjacentDisagreements`: count of benign-to-VUS, VUS-to-benign, VUS-to-oncogenic, and oncogenic-to-VUS disagreements in the collapsed 3x3 matrix
- `extremeDisagreements`: count of benign-to-oncogenic and oncogenic-to-benign disagreements in the collapsed 3x3 matrix
- `diagonalAgreement`: count of exact matches in the collapsed 3x3 matrix
- `exactScoreAgreement`: count of rows where the final numeric score exactly matched the reference score
- `exactScoreAgreementRate`: `exactScoreAgreement / numScoreTested`
- `scoreWithin2`: count of rows where the predicted overall score was within 2 points of the reference score
- `scoreWithin2Rate`: `scoreWithin2 / numScoreTested`
- `scoreGreaterThan2`: count of rows where the predicted overall score differed from the reference score by more than 2 points
- `scoreGreaterThan2Rate`: `scoreGreaterThan2 / numScoreTested`

### `evaluation/output/metrics-score-based.csv`

This file summarizes criteria and swimlane agreement.

Current contents include:

- input-row count
- unavailable-prediction count
- criteria-tested row count
- exact criteria-set match count and rate
- per-swimlane agreement rows for `population`, `computational`, `hotspots`, `predictive`, `om1`, `op2`, and `functional`

Per-swimlane rows report exact agreement, agreement within 2 points, larger deviations, and weighted agreement measures.

Rows without reference criteria are excluded from this file.

Field definitions:

- `numInputRows`: total number of input rows processed by the evaluator
- `numPredictionUnavailable`: number of rows where no final prediction was returned and the row was excluded from concordance metrics
- `numCriteriaTested`: number of rows with usable reference criteria and predicted criteria for score-based evaluation
- `exactCriteriaSetMatches`: number of rows where the full predicted criteria set exactly matched the reference criteria set
- `exactCriteriaSetAgreementRate`: `exactCriteriaSetMatches / numCriteriaTested`

Per-swimlane row definitions:

- `numTested`: number of rows included in comparison for that swimlane section
- `numCriteriaTested`: repeated total number of rows with usable criteria for convenience in each swimlane row
- `swimlane`: swimlane name being scored, currently one of `population`, `computational`, `hotspots`, `predictive`, `om1`, `op2`, or `functional`
- `same`: number of rows where the predicted swimlane score exactly matched the reference swimlane score
- `within2`: number of rows where the predicted swimlane score was within 2 points of the reference swimlane score
- `greaterThan2`: number of rows where the predicted swimlane score differed from the reference swimlane score by more than 2 points
- `sumCheck`: `same + within2 + greaterThan2`; used as a consistency check against `numTested`
- `exactAgreementRate`: `same / numTested`
- `acceptableAgreementRate`: `(same + within2) / numTested`
- `largeDeviationRate`: `greaterThan2 / numTested`
- `weightedScoreLinear`: `(same + 0.5 * within2) / numTested`
- `fuzzinessRatio`: `within2 / (same + within2)`; measures how much acceptable agreement depends on tolerance rather than exact matches
- `failureToExactRatio`: `greaterThan2 / same`; measures how often large failures occur relative to exact matches

## Unavailable Predictions

Rows where the final prediction is unavailable are still written to `oncogenicityPredictions.csv`, but they are excluded from both metrics files.

This keeps operational failures such as transient VEP outages from being miscounted as classification disagreements while still preserving them for manual review.