"""Shared gene-role classification helpers for evidence pipelines."""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from app.models.tumor_types import TumorType


GENE_ROLE_DATA_PATH = Path(__file__).resolve().parents[3] / "data" / "_Dict_Gene.csv"

GENE_ROLE_ONCOGENE = "oncogene"
GENE_ROLE_TSG = "tsg"
GENE_ROLE_BOTH = "oncogene_and_tsg"
GENE_ROLE_NEITHER = "neither"

DUAL_ROLE_TUMOR_TYPE_ROLE: dict[str, dict[str, str]] = {
    "GATA3": {
        "Peripheral T-Cell Lymphoma": GENE_ROLE_ONCOGENE,
        "T-Cell Acute Lymphoblastic Leukemia": GENE_ROLE_ONCOGENE,
        "Hodgkin Lymphoma": GENE_ROLE_ONCOGENE,
        "Neuroblastoma": GENE_ROLE_ONCOGENE,
        "T-Cell Lymphoblastic Lymphoma": GENE_ROLE_ONCOGENE,
        "Breast Cancer": GENE_ROLE_TSG,
        "Urothelial Carcinoma": GENE_ROLE_TSG,
        "Bladder Carcinoma": GENE_ROLE_TSG,
        "Renal Cell Carcinoma": GENE_ROLE_TSG,
        "Parathyroid Carcinoma": GENE_ROLE_TSG,
    }
}


@lru_cache(maxsize=1)
def load_gene_roles() -> dict[str, str | None]:
    roles: dict[str, str | None] = {}
    with GENE_ROLE_DATA_PATH.open(newline="", encoding="utf-8-sig") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            symbol = (row.get("NCBIgeneSymbol") or "").strip().upper()
            if not symbol:
                continue
            is_oncogene = (row.get("Oncogene") or "").strip().upper() == "YES"
            is_tsg = (row.get("TSG") or "").strip().upper() == "YES"
            if is_oncogene and is_tsg:
                roles[symbol] = GENE_ROLE_BOTH
            elif is_oncogene:
                roles[symbol] = GENE_ROLE_ONCOGENE
            elif is_tsg:
                roles[symbol] = GENE_ROLE_TSG
            else:
                roles[symbol] = None
    return roles


def classify_base_gene_role(gene: str | None) -> str | None:
    if not gene:
        return None
    role = load_gene_roles().get(gene.upper())
    if role is None:
        return GENE_ROLE_NEITHER
    return role


def resolve_gene_role(gene: str | None, tumor_type: TumorType | None) -> str | None:
    base_role = classify_base_gene_role(gene)
    if base_role in {GENE_ROLE_ONCOGENE, GENE_ROLE_TSG, GENE_ROLE_NEITHER, None}:
        return base_role
    if tumor_type is None:
        return None
    return DUAL_ROLE_TUMOR_TYPE_ROLE.get((gene or "").upper(), {}).get(tumor_type)
