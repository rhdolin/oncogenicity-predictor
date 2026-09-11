from fastapi.testclient import TestClient
import pytest

from app.main import app
from app.services.normalization import variant_normalizer


CLINGEN_SAMPLE_RECORD = {
    "@id": "http://reg.genome.network/allele/CA16602564",
    "genomicAlleles": [
        {
            "chromosome": "13",
            "coordinates": [
                {
                    "allele": "C",
                    "end": 28027222,
                    "referenceAllele": "A",
                    "start": 28027221,
                }
            ],
            "hgvs": ["NC_000013.11:g.28027222A>C", "CM000675.2:g.28027222A>C"],
            "referenceGenome": "GRCh38",
        },
        {
            "chromosome": "13",
            "coordinates": [
                {
                    "allele": "C",
                    "end": 28601359,
                    "referenceAllele": "A",
                    "start": 28601358,
                }
            ],
            "hgvs": ["NC_000013.10:g.28601359A>C", "CM000675.1:g.28601359A>C"],
            "referenceGenome": "GRCh37",
        },
    ],
    "transcriptAlleles": [
        {
            "geneSymbol": "FLT3",
            "geneNCBI_id": 2322,
            "hgvs": ["NM_004119.3:c.2073T>G", "NM_004119.2:c.2073T>G"],
            "genomeAlignments": [
                {
                    "referenceGenome": "GRCh37",
                    "chromosome": "13",
                    "start": 28601358,
                    "end": 28601359,
                }
            ],
            "proteinEffect": {
                "hgvs": "NP_004110.2:p.Phe691Leu",
                "hgvsWellDefined": "NP_004110.2:p.Phe691Leu",
            },
            "MANE": {
                "maneStatus": "MANE Select",
                "nucleotide": {
                    "RefSeq": {
                        "hgvs": "NM_004119.3:c.2073T>G",
                    }
                },
                "protein": {
                    "RefSeq": {
                        "hgvs": "NP_004110.2:p.Phe691Leu",
                    }
                },
            },
        }
    ],
}


client = TestClient(app)


@pytest.fixture(autouse=True)
def stub_clingen_lookup(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_fetch_allele_registry_record(submitted_variant: str) -> dict:
        assert submitted_variant
        return CLINGEN_SAMPLE_RECORD

    monkeypatch.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        fake_fetch_allele_registry_record,
    )


def test_root_exposes_docs() -> None:
    response = client.get("/")

    assert response.status_code == 200
    assert response.json()["docs"] == "/docs"


def test_predict_single_returns_normalized_variant() -> None:
    response = client.get(
        "/predictOncogenicity",
        params={"variant": "NM_004119.3:c.2073T>G"},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["submitted_variant"] == "NM_004119.3:c.2073T>G"
    assert body["identifiers"]["caid"] == "CA16602564"
    assert body["geneSymbol"] == "FLT3"
    assert body["geneNCBI_id"] == 2322
    assert body["genomic_hgvs"]["GRCh38"] == "NC_000013.11:g.28027222A>C"
    assert body["transcript_hgvs"]["canonical_b37"] == "NM_004119.3:c.2073T>G"
    assert body["transcript_hgvs"]["representative_transcript_hgvs"] == "NM_004119.3:c.2073T>G"
    assert body["protein"]["hgvs_1letter"] == "p.F691L"


def test_predict_batch_returns_normalized_variants() -> None:
    response = client.post(
        "/predictOncogenicity",
        json={"variants": ["NM_004119.3:c.2073T>G", "ENST00000241453.12:c.2073T>G"]},
    )

    assert response.status_code == 200
    body = response.json()
    assert len(body["normalized_variants"]) == 2
    assert body["normalized_variants"][1]["submitted_variant"] == "ENST00000241453.12:c.2073T>G"


def test_non_refseq_protein_effect_is_not_used() -> None:
    non_refseq_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                "geneSymbol": "FLT3",
                "geneNCBI_id": 2322,
                "hgvs": ["NM_004119.3:c.2073T>G"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "13",
                        "start": 28601358,
                        "end": 28601359,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "ENSP00000241453.7:p.Phe691Leu",
                    "hgvsWellDefined": "ENSP00000241453.7:p.Phe691Leu",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_004119.3:c.2073T>G",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "ENSP00000241453.7:p.Phe691Leu",
                        }
                    },
                },
            }
        ],
    }

    normalized_variant = variant_normalizer.normalize_variant("NM_004119.3:c.2073T>G")
    assert normalized_variant.protein.hgvs_protein_full == "NP_004110.2:p.Phe691Leu"

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: non_refseq_record
        normalized_variant = variant_normalizer.normalize_variant("NC_000009.12:g.130714395G>A")
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert normalized_variant.protein.np_accession is None
    assert normalized_variant.protein.hgvs_protein_full is None
    assert normalized_variant.protein.hgvs_3letter is None
    assert normalized_variant.protein.hgvs_1letter is None


def test_non_nc_nm_values_are_not_used() -> None:
    non_nc_nm_record = {
        **CLINGEN_SAMPLE_RECORD,
        "genomicAlleles": [
            {
                "chromosome": "13",
                "coordinates": [
                    {
                        "allele": "C",
                        "end": 28601359,
                        "referenceAllele": "A",
                        "start": 28601358,
                    }
                ],
                "hgvs": ["CM000675.1:g.28601359A>C"],
                "referenceGenome": "GRCh37",
            }
        ],
        "transcriptAlleles": [
            {
                "geneSymbol": "FLT3",
                "geneNCBI_id": 2322,
                "hgvs": ["NM_004119.3:c.2073T>G"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "13",
                        "start": 28601358,
                        "end": 28601359,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_004110.2:p.Phe691Leu",
                    "hgvsWellDefined": "NP_004110.2:p.Phe691Leu",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NR_130706.2:n.2139T>G",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_004110.2:p.Phe691Leu",
                        }
                    },
                },
            }
        ],
    }

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: non_nc_nm_record
        normalized_variant = variant_normalizer.normalize_variant("NC_000009.12:g.130714395G>A")
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert normalized_variant.genomic_hgvs.GRCh37 is None
    assert normalized_variant.coordinates.GRCh37.refseq is None
    assert normalized_variant.transcript_hgvs.mane_select_b38 is None
    assert normalized_variant.transcript_hgvs.canonical_b37 == "NM_004119.3:c.2073T>G"
    assert normalized_variant.transcript_hgvs.representative_transcript_hgvs == "NM_004119.3:c.2073T>G"


def test_representative_transcript_hgvs_falls_back_to_first_nm() -> None:
    transcript_fallback_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                "geneSymbol": "FLT3",
                "geneNCBI_id": 2322,
                "hgvs": ["NM_004119.2:c.2073T>G", "NR_130706.2:n.2139T>G"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh38",
                        "chromosome": "13",
                        "start": 28027221,
                        "end": 28027222,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_004110.2:p.Phe691Leu",
                    "hgvsWellDefined": "NP_004110.2:p.Phe691Leu",
                },
                "MANE": {
                    "maneStatus": "MANE Plus Clinical",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NR_130706.2:n.2139T>G",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_004110.2:p.Phe691Leu",
                        }
                    },
                },
            }
        ],
    }

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: transcript_fallback_record
        normalized_variant = variant_normalizer.normalize_variant("NC_000009.12:g.130714395G>A")
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert normalized_variant.transcript_hgvs.mane_select_b38 is None
    assert normalized_variant.transcript_hgvs.canonical_b37 is None
    assert normalized_variant.transcript_hgvs.representative_transcript_hgvs == "NM_004119.2:c.2073T>G"


def test_chr_x_and_chr_y_use_numeric_chrom_num() -> None:
    assert variant_normalizer._extract_coordinates(
        {
            "chromosome": "X",
            "coordinates": [{"allele": "A", "end": 32389644, "referenceAllele": "G"}],
        },
        "NC_000023.11:g.32389644G>A",
    ).chrom_num == "23"

    assert variant_normalizer._extract_coordinates(
        {
            "chromosome": "Y",
            "coordinates": [{"allele": "A", "end": 2787423, "referenceAllele": "G"}],
        },
        "NC_000024.10:g.2787423G>A",
    ).chrom_num == "24"

    mitochondrial_coordinates = variant_normalizer._extract_coordinates(
        {
            "chromosome": "MT",
            "coordinates": [{"allele": "G", "end": 3243, "referenceAllele": "A"}],
        },
        "NC_012920.1:g.3243A>G",
    )

    assert mitochondrial_coordinates.chrom == "chrM"
    assert mitochondrial_coordinates.chrom_num == "M"
