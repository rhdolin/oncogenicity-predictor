from fastapi.testclient import TestClient
import pytest

from app.main import app
from app.services.annotation import variant_annotator
from app.services.normalization import variant_normalizer


REAL_FETCH_VEP_ANNOTATION_RECORD = variant_annotator.fetch_vep_annotation_record


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

VEP_SAMPLE_RECORD = {
    "input": "NC_000013.11:g.28027222A>C",
    "allele_string": "A/C",
    "most_severe_consequence": "missense_variant",
    "transcript_consequences": [
        {
            "transcript_id": "NM_004119.3",
            "hgvsc": "NM_004119.3:c.2073T>G",
            "mane_select": "ENST00000241453.12",
            "mane": ["MANE_Select"],
            "consequence_terms": ["missense_variant"],
            "protein_start": 691,
            "protein_end": 691,
            "amino_acids": "F/L",
            "cadd_phred": 25.3,
            "cadd_raw": 4.12,
            "phylop100way_vertebrate": 7.89,
            "gene_symbol": "FLT3",
        },
        {
            "transcript_id": "NR_130706.2",
            "consequence_terms": ["non_coding_transcript_exon_variant"],
            "cadd_phred": 25.3,
            "cadd_raw": 4.12,
        },
    ],
    "colocated_variants": [
        {
            "frequencies": {
                "C": {
                    "gnomade_nfe": 0.0002,
                    "gnomade": 0.0001,
                    "gnomadg": 0.00015,
                }
            }
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
    monkeypatch.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: VEP_SAMPLE_RECORD,
    )


def test_root_exposes_docs() -> None:
    response = client.get("/")

    assert response.status_code == 200
    assert response.json()["docs"] == "/docs"


def test_predict_single_returns_annotated_variant() -> None:
    response = client.get(
        "/predictOncogenicity",
        params={"variant": "NM_004119.3:c.2073T>G"},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["normalizedVariant"]["submitted_variant"] == "NM_004119.3:c.2073T>G"
    assert body["normalizedVariant"]["identifiers"]["caid"] == "CA16602564"
    assert body["normalizedVariant"]["geneSymbol"] == "FLT3"
    assert body["normalizedVariant"]["geneNCBI_id"] == 2322
    assert body["normalizedVariant"]["genomic_hgvs"]["GRCh38"] == "NC_000013.11:g.28027222A>C"
    assert body["normalizedVariant"]["transcript_hgvs"]["mane_select_b38_source"] == "clingen"
    assert body["normalizedVariant"]["transcript_hgvs"]["canonical_b37"] == "NM_004119.3:c.2073T>G"
    assert body["normalizedVariant"]["transcript_hgvs"]["representative_transcript_hgvs"] == "NM_004119.3:c.2073T>G"
    assert body["normalizedVariant"]["protein"]["hgvs_1letter"] == "p.F691L"
    assert body["annotationStatus"] == "complete"
    assert body["annotationError"] is None
    assert body["basicAnnotation"]["mostSevereConsequence"] == "missense_variant"
    assert body["basicAnnotation"]["transcriptConsequences"] == [
        {
            "transcriptRefSeq": "NM_004119.3",
            "consequenceTerms": ["missense_variant"],
            "proteinStart": 691,
            "proteinEnd": 691,
            "aminoAcids": "F/L",
            "isManeSelect": True,
        }
    ]
    assert body["basicAnnotation"]["population"]["maxSubpopulationLabel"] == "gnomade_nfe"
    assert body["computationalAnnotation"]["cadd"]["phred"] == 25.3
    assert body["computationalAnnotation"]["phyloP100wayVertebrate"] == 7.89


def test_predict_batch_returns_annotated_variants() -> None:
    response = client.post(
        "/predictOncogenicity",
        json={"variants": ["NM_004119.3:c.2073T>G", "ENST00000241453.12:c.2073T>G"]},
    )

    assert response.status_code == 200
    body = response.json()
    assert len(body["annotated_variants"]) == 2
    assert body["annotated_variants"][1]["normalizedVariant"]["submitted_variant"] == "ENST00000241453.12:c.2073T>G"
    assert body["annotated_variants"][1]["annotationStatus"] == "complete"


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
    assert normalized_variant.coordinates.GRCh37 is not None
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


def test_vep_query_falls_back_from_grch38_to_mane_select() -> None:
    normalized_variant = variant_normalizer.normalize_variant("NM_004119.3:c.2073T>G")
    attempted_urls: list[str] = []

    class FakeResponse:
        def __init__(self, status_code: int, payload: list[dict]) -> None:
            self.status_code = status_code
            self._payload = payload

        def json(self) -> list[dict]:
            return self._payload

    def fake_httpx_get(url: str, **kwargs) -> FakeResponse:
        attempted_urls.append(url)
        if url.endswith("NC_000013.11%3Ag.28027222A%3EC"):
            return FakeResponse(404, [])
        if url.endswith("NM_004119.3%3Ac.2073T%3EG"):
            return FakeResponse(200, [VEP_SAMPLE_RECORD])
        return FakeResponse(500, [])

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(variant_annotator.httpx, "get", fake_httpx_get)
    try:
        record = REAL_FETCH_VEP_ANNOTATION_RECORD(normalized_variant)
    finally:
        monkeypatch_context.undo()

    assert record == VEP_SAMPLE_RECORD
    assert attempted_urls == [
        "https://rest.ensembl.org/vep/human/hgvs/NC_000013.11%3Ag.28027222A%3EC",
        "https://rest.ensembl.org/vep/human/hgvs/NM_004119.3%3Ac.2073T%3EG",
    ]


def test_annotate_variant_backfills_mane_from_vep_when_clingen_lacks_it() -> None:
    clingen_without_mane_record = {
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

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        lambda submitted_variant: clingen_without_mane_record,
    )
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: VEP_SAMPLE_RECORD,
    )
    try:
        annotated_variant = variant_annotator.annotate_variant(
            variant_normalizer.normalize_variant("NC_000013.11:g.28027222A>C")
        )
    finally:
        monkeypatch_context.undo()

    assert annotated_variant.normalizedVariant.transcript_hgvs.mane_select_b38 == "NM_004119.3:c.2073T>G"
    assert annotated_variant.normalizedVariant.transcript_hgvs.mane_select_b38_source == "vep"


def test_predict_single_returns_failed_annotation_payload_when_vep_fails() -> None:
    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: (_ for _ in ()).throw(
            variant_annotator.VariantAnnotationError(
                "VEP annotation failed for all supported query forms.",
                [
                    "NC_000013.11:g.28027222A>C",
                    "NM_004119.3:c.2073T>G",
                ],
            )
        ),
    )
    try:
        response = client.get(
            "/predictOncogenicity",
            params={"variant": "NM_004119.3:c.2073T>G"},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    assert body["normalizedVariant"]["submitted_variant"] == "NM_004119.3:c.2073T>G"
    assert body["annotationStatus"] == "failed"
    assert body["annotationError"] == {
        "source": "vep",
        "message": "VEP annotation failed for all supported query forms.",
        "attemptedQueries": [
            "NC_000013.11:g.28027222A>C",
            "NM_004119.3:c.2073T>G",
        ],
    }
    assert body["basicAnnotation"] is None
    assert body["computationalAnnotation"] is None


def test_predict_batch_returns_mixed_success_and_failed_annotation_payloads() -> None:
    def fake_fetch_vep_annotation_record(normalized_variant):
        if normalized_variant.submitted_variant == "ENST00000241453.12:c.2073T>G":
            raise variant_annotator.VariantAnnotationError(
                "VEP annotation failed for all supported query forms.",
                [
                    "NC_000013.11:g.28027222A>C",
                    "NM_004119.3:c.2073T>G",
                ],
            )
        return VEP_SAMPLE_RECORD

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        fake_fetch_vep_annotation_record,
    )
    try:
        response = client.post(
            "/predictOncogenicity",
            json={"variants": ["NM_004119.3:c.2073T>G", "ENST00000241453.12:c.2073T>G"]},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    assert body["annotated_variants"][0]["annotationStatus"] == "complete"
    assert body["annotated_variants"][1]["annotationStatus"] == "failed"
    assert body["annotated_variants"][1]["basicAnnotation"] is None
    assert body["annotated_variants"][1]["computationalAnnotation"] is None


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
