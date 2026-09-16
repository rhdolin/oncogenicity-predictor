"""ClinMAVE-backed functional evidence scoring."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from app.models.annotated_variant import AnnotatedVariant
from app.models.prediction import EvidenceResult
from app.models.tumor_types import TumorType
from app.services.evidence.gene_roles import (
    GENE_ROLE_BOTH,
    GENE_ROLE_ONCOGENE,
    GENE_ROLE_TSG,
    classify_base_gene_role,
    resolve_gene_role,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
CLINMAVE_DATA_DIR = REPO_ROOT / "data" / "clinmave"

CLASSIFICATION_NORMAL = "Functionally normal"
CLASSIFICATION_GOF = "Gain-of-function"
CLASSIFICATION_LOF = "Loss-of-function"

IDENTIFIER_TRANSCRIPT_HGVS_RE = re.compile(r"^(NM_[0-9]+\.[0-9]+)\([^)]+\):(c\.[^ ]+)")


@dataclass(slots=True)
class ClinMaveRow:
    identifier: str
    transcript_hgvs: str
    gene: str
    functional_classification: str
    functional_description: str
    molecular_consequence: str
    phenotype: str
    dataset_id: str
    score: str
    publication: str


_GENE_INDEX_CACHE: dict[str, dict[str, list[ClinMaveRow]]] = {}


def build_functional_evidence(
    annotated_variant: AnnotatedVariant,
    tumor_type: TumorType | None = None,
) -> EvidenceResult:
    if annotated_variant.annotationStatus != "complete":
        return _build_not_available_result(
            statement="Functional evidence could not be evaluated because annotation data was unavailable.",
            matched_data=None,
            absent_reason="error",
        )

    normalized_variant = annotated_variant.normalizedVariant
    gene = (normalized_variant.geneSymbol or "").upper()
    mane_select_b38 = normalized_variant.transcript_hgvs.mane_select_b38

    if not gene:
        return _build_not_available_result(
            statement="Functional evidence not available: no gene symbol was resolved for this variant.",
            matched_data={"tumorType": tumor_type},
        )

    if not mane_select_b38:
        return _build_not_available_result(
            statement=(
                "Functional evidence not available: no MANE transcript HGVS "
                f"representation was resolved for gene {gene}."
            ),
            matched_data={"gene": gene, "tumorType": tumor_type},
        )

    gene_index = _get_gene_index(gene)
    if gene_index is None:
        return _build_not_available_result(
            statement=(
                "Functional evidence not available: gene "
                f"{gene} is not included in the retained ClinMAVE dataset."
            ),
            matched_data={
                "gene": gene,
                "tumorType": tumor_type,
                "maneSelectB38": mane_select_b38,
            },
        )

    matched_rows = gene_index.get(mane_select_b38)
    if matched_rows is None:
        return _build_not_available_result(
            statement=(
                "Functional evidence not available: this variant was not found "
                f"in the retained ClinMAVE dataset for gene {gene}."
            ),
            matched_data={
                "gene": gene,
                "tumorType": tumor_type,
                "maneSelectB38": mane_select_b38,
            },
        )

    classifications = {
        row.functional_classification
        for row in matched_rows
        if row.functional_classification
    }
    if len(classifications) != 1:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "ClinMAVE contains conflicting functional classifications for "
                f"gene {gene} and variant {mane_select_b38}; no functional "
                "evidence rule is applied."
            ),
            status="applied",
            source="clinmave",
            matchedData=_build_matched_data(matched_rows, gene, tumor_type, None),
        )

    matched_row = matched_rows[-1]
    resolved_role = resolve_gene_role(gene, tumor_type)
    if resolved_role is None:
        if classify_base_gene_role(gene) == GENE_ROLE_BOTH:
            if tumor_type is None:
                statement = (
                    "Functional evidence not available: gene "
                    f"{gene} requires tumor type to resolve functional interpretation."
                )
            else:
                statement = (
                    "Functional evidence not available: no ClinMAVE role "
                    f"mapping is defined for gene {gene} in tumor type {tumor_type}."
                )
        else:
            statement = (
                "Functional evidence not available: no supported gene role "
                f"mapping is defined for gene {gene}."
            )
        return _build_not_available_result(
            statement=statement,
            matched_data=_build_matched_data(
                matched_rows,
                gene,
                tumor_type,
                resolved_role,
            ),
        )

    classification = matched_row.functional_classification
    matched_data = _build_matched_data(matched_rows, gene, tumor_type, resolved_role)

    if resolved_role == GENE_ROLE_ONCOGENE and classification == CLASSIFICATION_GOF:
        return EvidenceResult(
            score=4,
            evidenceCode="OS2",
            evidenceStatement=(
                f"ClinMAVE shows gain-of-function for this variant in oncogene {gene}."
            ),
            status="applied",
            source="clinmave",
            matchedData=matched_data,
        )

    if resolved_role == GENE_ROLE_TSG and classification == CLASSIFICATION_LOF:
        return EvidenceResult(
            score=4,
            evidenceCode="OS2",
            evidenceStatement=(
                f"ClinMAVE shows loss-of-function for this variant in tumor suppressor gene {gene}."
            ),
            status="applied",
            source="clinmave",
            matchedData=matched_data,
        )

    if resolved_role == GENE_ROLE_ONCOGENE and classification == CLASSIFICATION_NORMAL:
        return EvidenceResult(
            score=-4,
            evidenceCode="SBS2",
            evidenceStatement=(
                f"ClinMAVE shows functionally normal activity for this variant in oncogene {gene}."
            ),
            status="applied",
            source="clinmave",
            matchedData=matched_data,
        )

    if resolved_role == GENE_ROLE_TSG and classification == CLASSIFICATION_NORMAL:
        return EvidenceResult(
            score=-4,
            evidenceCode="SBS2",
            evidenceStatement=(
                "ClinMAVE shows functionally normal activity for this "
                f"variant in tumor suppressor gene {gene}."
            ),
            status="applied",
            source="clinmave",
            matchedData=matched_data,
        )

    if resolved_role == GENE_ROLE_ONCOGENE and classification == CLASSIFICATION_LOF:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "ClinMAVE shows loss-of-function for this variant in oncogene "
                f"{gene}; no functional evidence rule is applied."
            ),
            status="applied",
            source="clinmave",
            matchedData=matched_data,
        )

    if resolved_role == GENE_ROLE_TSG and classification == CLASSIFICATION_GOF:
        return EvidenceResult(
            score=0,
            evidenceCode=None,
            evidenceStatement=(
                "ClinMAVE shows gain-of-function for this variant in tumor "
                f"suppressor gene {gene}; no functional evidence rule is applied."
            ),
            status="applied",
            source="clinmave",
            matchedData=matched_data,
        )

    return _build_not_available_result(
        statement=(
            "Functional evidence not available: ClinMAVE classification "
            f"{classification!r} is not supported for gene {gene}."
        ),
        matched_data=matched_data,
    )


def _build_not_available_result(
    statement: str,
    matched_data: dict[str, str | None] | None,
    absent_reason: str = "unsupported",
) -> EvidenceResult:
    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement=statement,
        status="not_available",
        source="clinmave",
        matchedData=matched_data,
        dataAbsentReason=absent_reason,
    )


def _build_matched_data(
    matched_rows: list[ClinMaveRow],
    gene: str,
    tumor_type: TumorType | None,
    resolved_role: str | None,
) -> dict[str, str | None]:
    matched_row = matched_rows[-1]
    return {
        "gene": gene,
        "tumorType": tumor_type,
        "resolvedGeneRole": resolved_role,
        "clinmaveIdentifier": matched_row.identifier,
        "maneSelectB38": matched_row.transcript_hgvs,
        "functionalClassification": matched_row.functional_classification,
        "functionalDescription": matched_row.functional_description,
        "molecularConsequence": matched_row.molecular_consequence,
        "phenotype": matched_row.phenotype,
        "datasetId": matched_row.dataset_id,
        "score": matched_row.score,
        "publication": matched_row.publication,
        "matchingRowCount": len(matched_rows),
        "matchingDatasetIds": [row.dataset_id for row in matched_rows],
        "matchingFunctionalClassifications": [
            row.functional_classification for row in matched_rows
        ],
    }


def _get_gene_index(gene: str) -> dict[str, list[ClinMaveRow]] | None:
    gene = gene.upper()
    cached = _GENE_INDEX_CACHE.get(gene)
    if cached is not None:
        return cached

    path = CLINMAVE_DATA_DIR / f"variants.{gene}.csv"
    if not path.exists():
        return None

    rows_by_transcript_hgvs: dict[str, list[ClinMaveRow]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        identifier_column = (reader.fieldnames or ["Identifier"])[0]
        for row in reader:
            identifier = (row.get(identifier_column) or "").strip()
            transcript_hgvs = _normalize_identifier_to_transcript_hgvs(identifier)
            if transcript_hgvs is None:
                continue
            rows_by_transcript_hgvs.setdefault(transcript_hgvs, []).append(
                ClinMaveRow(
                    identifier=identifier,
                    transcript_hgvs=transcript_hgvs,
                    gene=(row.get("Gene name") or gene).strip().upper(),
                    functional_classification=(
                        row.get("Functional classification") or ""
                    ).strip(),
                    functional_description=(
                        row.get("Functional description") or ""
                    ).strip(),
                    molecular_consequence=(
                        row.get("Molecular consequence") or ""
                    ).strip(),
                    phenotype=(row.get("Phenotype") or "").strip(),
                    dataset_id=(row.get("Dataset ID") or "").strip(),
                    score=(row.get("Score") or "").strip(),
                    publication=(row.get("Publication") or "").strip(),
                )
            )

    _GENE_INDEX_CACHE[gene] = rows_by_transcript_hgvs
    return rows_by_transcript_hgvs


def _normalize_identifier_to_transcript_hgvs(identifier: str) -> str | None:
    match = IDENTIFIER_TRANSCRIPT_HGVS_RE.match(identifier)
    if match is None:
        return None
    return f"{match.group(1)}:{match.group(2)}"
