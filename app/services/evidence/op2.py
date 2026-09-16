"""OP2 evidence scoring for curated single-etiology tumor contexts."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from app.models.annotated_variant import AnnotatedVariant, TranscriptConsequence
from app.models.prediction import EvidenceResult
from app.models.tumor_types import TumorType


REPO_ROOT = Path(__file__).resolve().parents[3]
OP2_RULES_PATH = REPO_ROOT / "data" / "op2_rules.csv"


class Op2RulesError(Exception):
    """Raised when the local OP2 rules table cannot be parsed."""


@dataclass(frozen=True, slots=True)
class Op2Rule:
    tumor_type: str
    gene_symbol: str
    mane_protein_hgvs: str | None


def build_op2_evidence(
    annotated_variant: AnnotatedVariant,
    tumor_type: TumorType | None = None,
) -> EvidenceResult:
    """Evaluate OP2 using the curated local tumor-type rules table."""
    if annotated_variant.annotationStatus != "complete" or annotated_variant.basicAnnotation is None:
        return _build_not_available_result(
            statement="OP2 evidence could not be evaluated because annotation data was unavailable.",
            matched_data=None,
            absent_reason="error",
        )

    gene_symbol = (annotated_variant.normalizedVariant.geneSymbol or "").upper()
    matched_data = {
        "tumorType": tumor_type,
        "gene": gene_symbol or None,
    }

    if tumor_type is None:
        return _build_not_available_result(
            statement="OP2 evidence not available: tumor type is required for this rule.",
            matched_data=matched_data,
        )

    if not gene_symbol:
        return _build_not_available_result(
            statement="OP2 evidence could not be evaluated because the normalized gene symbol was unavailable.",
            matched_data=matched_data,
            absent_reason="unknown",
        )

    try:
        rules = _load_op2_rules()
    except (FileNotFoundError, Op2RulesError) as exc:
        return _build_not_available_result(
            statement=(
                "OP2 evidence could not be evaluated because the local OP2 rules "
                f"table was unavailable or unreadable: {exc}"
            ),
            matched_data=matched_data,
            absent_reason="error",
        )

    candidate_rules = [
        rule
        for rule in rules
        if rule.tumor_type == tumor_type and rule.gene_symbol == gene_symbol
    ]
    matched_data["candidateRuleCount"] = len(candidate_rules)

    if not candidate_rules:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement="OP2 evidence did not meet current scoring criteria.",
            status="applied",
            source="op2",
            matchedData=matched_data,
            dataAbsentReason=None,
        )

    specific_rules = [rule for rule in candidate_rules if rule.mane_protein_hgvs]
    broad_rules = [rule for rule in candidate_rules if not rule.mane_protein_hgvs]

    if specific_rules:
        mane_protein_hgvs = _extract_mane_protein_hgvs(
            annotated_variant.basicAnnotation.transcriptConsequences
        )
        matched_data["maneProteinHgvs"] = mane_protein_hgvs
        matched_data["specificRuleProteinHgvs"] = [
            rule.mane_protein_hgvs for rule in specific_rules if rule.mane_protein_hgvs
        ]

        if mane_protein_hgvs is None:
            return _build_not_available_result(
                statement=(
                    "OP2 evidence not available: a MANE Select protein consequence was "
                    f"required for tumor type {tumor_type} and gene {gene_symbol}."
                ),
                matched_data=matched_data,
            )

        for rule in specific_rules:
            if rule.mane_protein_hgvs == mane_protein_hgvs:
                return EvidenceResult(
                    score=1,
                    evidenceCode="OP2",
                    evidenceStatement=(
                        "Somatic variant in a gene associated with a malignancy "
                        "with a curated single genetic etiology context."
                    ),
                    status="applied",
                    source="op2",
                    matchedData={
                        **matched_data,
                        "matchedRule": {
                            "tumorType": rule.tumor_type,
                            "geneSymbol": rule.gene_symbol,
                            "maneProteinHgvs": rule.mane_protein_hgvs,
                        },
                    },
                    dataAbsentReason=None,
                )

    if broad_rules:
        matched_rule = broad_rules[0]
        return EvidenceResult(
            score=1,
            evidenceCode="OP2",
            evidenceStatement=(
                "Somatic variant in a gene associated with a malignancy with a curated "
                "single genetic etiology context."
            ),
            status="applied",
            source="op2",
            matchedData={
                **matched_data,
                "matchedRule": {
                    "tumorType": matched_rule.tumor_type,
                    "geneSymbol": matched_rule.gene_symbol,
                    "maneProteinHgvs": matched_rule.mane_protein_hgvs,
                },
            },
            dataAbsentReason=None,
        )

    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement="OP2 evidence did not meet current scoring criteria.",
        status="applied",
        source="op2",
        matchedData=matched_data,
        dataAbsentReason=None,
    )


def _build_not_available_result(
    statement: str,
    matched_data: dict | None,
    absent_reason: str = "unsupported",
) -> EvidenceResult:
    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement=statement,
        status="not_available",
        source="op2",
        matchedData=matched_data,
        dataAbsentReason=absent_reason,
    )


def _extract_mane_protein_hgvs(
    transcript_consequences: list[TranscriptConsequence],
) -> str | None:
    for consequence in transcript_consequences:
        if consequence.isManeSelect:
            return consequence.proteinHgvs
    return None


@lru_cache(maxsize=1)
def _load_op2_rules() -> tuple[Op2Rule, ...]:
    if not OP2_RULES_PATH.exists():
        raise FileNotFoundError(OP2_RULES_PATH)

    rules: list[Op2Rule] = []
    with OP2_RULES_PATH.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        expected_columns = {"tumorType", "geneSymbol", "maneProteinHgvs"}
        if set(reader.fieldnames or []) != expected_columns:
            raise Op2RulesError(
                "expected columns tumorType,geneSymbol,maneProteinHgvs"
            )

        for row in reader:
            tumor_type = (row.get("tumorType") or "").strip()
            gene_symbol = (row.get("geneSymbol") or "").strip().upper()
            mane_protein_hgvs = (row.get("maneProteinHgvs") or "").strip() or None
            if not tumor_type or not gene_symbol:
                continue
            rules.append(
                Op2Rule(
                    tumor_type=tumor_type,
                    gene_symbol=gene_symbol,
                    mane_protein_hgvs=mane_protein_hgvs,
                )
            )

    return tuple(rules)
