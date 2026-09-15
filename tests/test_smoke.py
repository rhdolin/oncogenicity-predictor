from fastapi.testclient import TestClient
import pytest
from datetime import datetime

from app.main import app
from app.models.annotated_variant import AnnotatedVariant, BasicAnnotation, TranscriptConsequence
from app.models.normalized_variant import NormalizationMetadata, NormalizedVariant
from app.services.annotation import variant_annotator
from app.services.evidence import functional
from app.services.evidence import gene_roles
from app.services.evidence import hotspots
from app.services.evidence import predictive
from app.services.normalization import variant_normalizer


REAL_FETCH_VEP_ANNOTATION_RECORD = variant_annotator.fetch_vep_annotation_record


FAKE_HOTSPOT_INDEX = hotspots.HotspotIndex(
    snv_by_gene={
        "FLT3": (
            hotspots.SnpHotspotRecord(
                gene="FLT3",
                position=691,
                ref_amino_acid="F",
                alt_amino_acid="L",
                mutation_count=75,
                variant_count=12,
            ),
        )
    },
    indel_by_gene={},
)

EMPTY_HOTSPOT_INDEX = hotspots.HotspotIndex(snv_by_gene={}, indel_by_gene={})


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
            "hgvsp": "NP_004110.2:p.Phe691Leu",
            "mane_select": "ENST00000241453.12",
            "mane": ["MANE_Select"],
            "consequence_terms": ["missense_variant"],
            "protein_start": 691,
            "protein_end": 691,
            "amino_acids": "F/L",
            "cadd_phred": 25.3,
            "cadd_raw": 4.12,
            "fathmm-xf_coding_pred": "D",
            "fathmm-xf_coding_score": 0.88,
            "fathmm-xf_coding_rankscore": 0.91,
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


def _component_by_code(body: dict, code: str) -> dict:
    return next(
        component
        for component in body["component"]
        if component["code"]["coding"][0]["code"] == code
    )


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
    monkeypatch.setattr(hotspots, "_get_hotspot_index", lambda: EMPTY_HOTSPOT_INDEX)
    monkeypatch.setattr(predictive, "search_clinvar_variation_ids", lambda query, retmax=100: [])
    monkeypatch.setattr(predictive, "fetch_clinvar_summaries", lambda variation_ids: {})
    functional._GENE_INDEX_CACHE.clear()
    gene_roles.load_gene_roles.cache_clear()


def test_root_exposes_docs() -> None:
    response = client.get("/")

    assert response.status_code == 200
    assert response.json()["docs"] == "/docs"


def test_predict_single_returns_observation() -> None:
    response = client.get(
        "/predictOncogenicity",
        params={"variant": "NM_004119.3:c.2073T>G"},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["resourceType"] == "Observation"
    assert body["status"] == "final"
    assert body["issued"].endswith("Z")
    datetime.fromisoformat(body["issued"].replace("Z", "+00:00"))
    assert body["code"]["coding"][0]["code"] == "oncogenicity-prediction"
    assert body["valueInteger"] == 2
    assert body["interpretation"] == []
    assert body["extension"] == [
        {
            "url": "https://oncogenicity-predictor.example/fhir/StructureDefinition/submitted-variant",
            "valueString": "NM_004119.3:c.2073T>G",
        }
    ]

    population = _component_by_code(body, "population-evidence")
    assert population["valueInteger"] == 1
    assert population["interpretation"] == [
        {
            "coding": [
                {
                    "system": "https://oncogenicity-predictor.example/fhir/CodeSystem/temp-codes",
                    "code": "OP4",
                    "display": "OP4",
                }
            ],
            "text": "Present at low frequency in gnomAD (≤1%; observed 0.02%).",
        }
    ]
    assert "dataAbsentReason" not in population

    computational = _component_by_code(body, "computational-evidence")
    assert computational["valueInteger"] == 1
    assert computational["interpretation"][0]["coding"][0]["code"] == "OP1"
    assert computational["interpretation"][0]["text"] == (
        "CADD supports oncogenicity for this variant (PHRED 25.3; most severe consequence missense_variant)."
    )
    assert "dataAbsentReason" not in computational

    hotspots = _component_by_code(body, "hotspots-evidence")
    assert hotspots["valueInteger"] == 0
    assert hotspots["interpretation"] == [
        {
            "coding": [],
            "text": "Hotspots evidence did not meet current scoring criteria.",
        }
    ]
    assert "dataAbsentReason" not in hotspots

    predictive = _component_by_code(body, "predictive-evidence")
    assert predictive["valueInteger"] == 0
    assert predictive["interpretation"][0]["text"] == (
        "Predictive evidence did not meet current scoring criteria."
    )
    assert "dataAbsentReason" not in predictive

    functional = _component_by_code(body, "functional-evidence")
    assert functional["dataAbsentReason"]["coding"][0]["code"] == "unsupported"


def test_annotate_single_returns_annotated_variant() -> None:
    response = client.get(
        "/annotateVariant",
        params={"variant": "NM_004119.3:c.2073T>G"},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["normalizedVariant"]["submitted_variant"] == "NM_004119.3:c.2073T>G"
    assert body["annotationStatus"] == "complete"
    assert body["basicAnnotation"]["mostSevereConsequence"] == "missense_variant"
    transcript_consequence = body["basicAnnotation"]["transcriptConsequences"][0]
    assert transcript_consequence["proteinHgvs"] == "p.F691L"
    assert transcript_consequence["proteinEventType"] == "substitution"
    assert transcript_consequence["rawProteinHgvs"] == "NP_004110.2:p.Phe691Leu"
    assert body["computationalAnnotation"]["fathmmXfCoding"]["prediction"] == "D"


def test_summarize_evidence_returns_raw_summary() -> None:
    response = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_004119.3:c.2073T>G"},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["overallScore"] == 2
    assert "overallClassification" in body
    assert body["overallClassification"] is None
    assert body["oncogenicityEvidence"]["population"]["evidenceCode"] == "OP4"
    assert body["oncogenicityEvidence"]["population"]["matchedData"]["effectiveAf"] == 0.0002
    assert body["oncogenicityEvidence"]["computational"]["evidenceCode"] == "OP1"
    assert body["oncogenicityEvidence"]["computational"]["matchedData"]["fathmmXfCodingPrediction"] == "D"
    assert body["oncogenicityEvidence"]["hotspots"]["evidenceCode"] is None
    assert body["oncogenicityEvidence"]["hotspots"]["matchedData"]["gene"] == "FLT3"
    assert body["oncogenicityEvidence"]["predictive"]["evidenceCode"] is None
    assert body["oncogenicityEvidence"]["predictive"]["matchedData"]["proteinChangeOneLetter"] == "F691L"


def test_functional_gene_not_in_retained_panel_is_not_available() -> None:
    response = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_004119.3:c.2073T>G"},
    )

    assert response.status_code == 200
    body = response.json()
    functional_evidence = body["oncogenicityEvidence"]["functional"]
    assert functional_evidence["status"] == "not_available"
    assert functional_evidence["matchedData"]["gene"] == "FLT3"


def test_functional_atm_match_applies_os2(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    clinmave_dir = tmp_path / "clinmave"
    clinmave_dir.mkdir()
    (clinmave_dir / "variants.ATM.csv").write_text(
        '"Identifier","Chrom","Position","Ref/Alt","Gene name","Score","Dataset ID","Molecular consequence","Functional description","Phenotype","Cross-assay hits","Functional classification","ClinVar information","Population frequency","TCGA summary","MAVE technique","Mutagenesis strategy","Publication"\n'  # noqa: E501
        '"NM_000051.4(ATM):c.283C>T (p.Gln95Ter)","chr11","108229275","C/T","ATM","-0.315","dataset0078","Nonsense","Reduced ATM function in model system","ATM-mediated anti-cell growth with Olaparib","5","Loss-of-function","Pathogenic/Likely pathogenic","6.19792e-07","NA","CRISPR-Based Genome Editing","Base editing","33606978"\n',  # noqa: E501
        encoding="utf-8",
    )

    monkeypatch.setattr(functional, "CLINMAVE_DATA_DIR", clinmave_dir)
    functional._GENE_INDEX_CACHE.clear()

    atm_record = {
        "@id": "http://reg.genome.network/allele/CA164660",
        "genomicAlleles": [
            {
                "chromosome": "11",
                "coordinates": [
                    {
                        "allele": "T",
                        "end": 108229275,
                        "referenceAllele": "C",
                        "start": 108229274,
                    }
                ],
                "hgvs": ["NC_000011.10:g.108229275C>T"],
                "referenceGenome": "GRCh38",
            }
        ],
        "transcriptAlleles": [
            {
                "geneSymbol": "ATM",
                "geneNCBI_id": 472,
                "hgvs": ["NM_000051.4:c.283C>T"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "11",
                        "start": 108100001,
                        "end": 108100002,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_000042.3:p.Gln95Ter",
                    "hgvsWellDefined": "NP_000042.3:p.Gln95Ter",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_000051.4:c.283C>T",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_000042.3:p.Gln95Ter",
                        }
                    },
                },
            }
        ],
    }

    atm_vep_record = {
        **VEP_SAMPLE_RECORD,
        "input": "NC_000011.10:g.108229275C>T",
        "most_severe_consequence": "stop_gained",
        "transcript_consequences": [
            {
                "transcript_id": "NM_000051.4",
                "hgvsc": "NM_000051.4:c.283C>T",
                "hgvsp": "NP_000042.3:p.Gln95Ter",
                "mane_select": "ENST00000278616.10",
                "mane": ["MANE_Select"],
                "consequence_terms": ["stop_gained"],
                "protein_start": 95,
                "protein_end": 95,
                "amino_acids": "Q/*",
                "cadd_phred": 36.0,
                "cadd_raw": 7.62,
                "fathmm-xf_coding_pred": "N",
                "fathmm-xf_coding_score": 0.09,
                "fathmm-xf_coding_rankscore": 0.18,
                "phylop100way_vertebrate": 5.8,
                "gene_symbol": "ATM",
            }
        ],
        "colocated_variants": [],
    }

    monkeypatch.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        lambda submitted_variant: atm_record,
    )
    monkeypatch.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: atm_vep_record,
    )

    response = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_000051.4:c.283C>T"},
    )

    assert response.status_code == 200
    body = response.json()
    functional_evidence = body["oncogenicityEvidence"]["functional"]
    assert functional_evidence["status"] == "applied"
    assert functional_evidence["score"] == 4
    assert functional_evidence["evidenceCode"] == "OS2"
    assert functional_evidence["matchedData"]["maneSelectB38"] == "NM_000051.4:c.283C>T"


def test_functional_atm_normal_match_applies_sbs2(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    clinmave_dir = tmp_path / "clinmave"
    clinmave_dir.mkdir()
    (clinmave_dir / "variants.ATM.csv").write_text(
        '"Identifier","Chrom","Position","Ref/Alt","Gene name","Score","Dataset ID","Molecular consequence","Functional description","Phenotype","Cross-assay hits","Functional classification","ClinVar information","Population frequency","TCGA summary","MAVE technique","Mutagenesis strategy","Publication"\n'  # noqa: E501
        '"NM_000051.4(ATM):c.283C>T (p.Gln95Ter)","chr11","108229275","C/T","ATM","-0.315","dataset0078","Nonsense","Neutral ATM function in model system","ATM-mediated anti-cell growth with Olaparib","5","Functionally normal","Pathogenic/Likely pathogenic","6.19792e-07","NA","CRISPR-Based Genome Editing","Base editing","33606978"\n',  # noqa: E501
        encoding="utf-8",
    )

    monkeypatch.setattr(functional, "CLINMAVE_DATA_DIR", clinmave_dir)
    functional._GENE_INDEX_CACHE.clear()

    atm_record = {
        "@id": "http://reg.genome.network/allele/CA164660",
        "genomicAlleles": [
            {
                "chromosome": "11",
                "coordinates": [
                    {
                        "allele": "T",
                        "end": 108229275,
                        "referenceAllele": "C",
                        "start": 108229274,
                    }
                ],
                "hgvs": ["NC_000011.10:g.108229275C>T"],
                "referenceGenome": "GRCh38",
            }
        ],
        "transcriptAlleles": [
            {
                "geneSymbol": "ATM",
                "geneNCBI_id": 472,
                "hgvs": ["NM_000051.4:c.283C>T"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "11",
                        "start": 108100001,
                        "end": 108100002,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_000042.3:p.Gln95Ter",
                    "hgvsWellDefined": "NP_000042.3:p.Gln95Ter",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_000051.4:c.283C>T",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_000042.3:p.Gln95Ter",
                        }
                    },
                },
            }
        ],
    }

    atm_vep_record = {
        **VEP_SAMPLE_RECORD,
        "input": "NC_000011.10:g.108229275C>T",
        "most_severe_consequence": "stop_gained",
        "transcript_consequences": [
            {
                "transcript_id": "NM_000051.4",
                "hgvsc": "NM_000051.4:c.283C>T",
                "hgvsp": "NP_000042.3:p.Gln95Ter",
                "mane_select": "ENST00000278616.10",
                "mane": ["MANE_Select"],
                "consequence_terms": ["stop_gained"],
                "protein_start": 95,
                "protein_end": 95,
                "amino_acids": "Q/*",
                "cadd_phred": 36.0,
                "cadd_raw": 7.62,
                "fathmm-xf_coding_pred": "N",
                "fathmm-xf_coding_score": 0.09,
                "fathmm-xf_coding_rankscore": 0.18,
                "phylop100way_vertebrate": 5.8,
                "gene_symbol": "ATM",
            }
        ],
        "colocated_variants": [],
    }

    monkeypatch.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        lambda submitted_variant: atm_record,
    )
    monkeypatch.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: atm_vep_record,
    )

    response = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_000051.4:c.283C>T"},
    )

    assert response.status_code == 200
    body = response.json()
    functional_evidence = body["oncogenicityEvidence"]["functional"]
    assert functional_evidence["status"] == "applied"
    assert functional_evidence["score"] == -4
    assert functional_evidence["evidenceCode"] == "SBS2"


def test_functional_gata3_requires_tumor_type(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    clinmave_dir = tmp_path / "clinmave"
    clinmave_dir.mkdir()
    (clinmave_dir / "variants.GATA3.csv").write_text(
        '"Identifier","Chrom","Position","Ref/Alt","Gene name","Score","Dataset ID","Molecular consequence","Functional description","Phenotype","Cross-assay hits","Functional classification","ClinVar information","Population frequency","TCGA summary","MAVE technique","Mutagenesis strategy","Publication"\n'  # noqa: E501
        '"NM_001002295.2(GATA3):c.1A>G (p.Met1Val)","chr10","8045419","A/G","GATA3","1.0","dataset0001","Missense","Activating effect in model system","example phenotype","1","Gain-of-function","","","NA","Deep Mutational Scanning","Saturation mutagenesis","12345678"\n',  # noqa: E501
        encoding="utf-8",
    )

    monkeypatch.setattr(functional, "CLINMAVE_DATA_DIR", clinmave_dir)
    functional._GENE_INDEX_CACHE.clear()

    gata3_record = {
        "@id": "http://reg.genome.network/allele/CATEST",
        "genomicAlleles": [
            {
                "chromosome": "10",
                "coordinates": [
                    {
                        "allele": "G",
                        "end": 8045419,
                        "referenceAllele": "A",
                        "start": 8045418,
                    }
                ],
                "hgvs": ["NC_000010.11:g.8045419A>G"],
                "referenceGenome": "GRCh38",
            }
        ],
        "transcriptAlleles": [
            {
                "geneSymbol": "GATA3",
                "geneNCBI_id": 2625,
                "hgvs": ["NM_001002295.2:c.1A>G"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "10",
                        "start": 8116649,
                        "end": 8116650,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_001002295.1:p.Met1Val",
                    "hgvsWellDefined": "NP_001002295.1:p.Met1Val",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_001002295.2:c.1A>G",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_001002295.1:p.Met1Val",
                        }
                    },
                },
            }
        ],
    }

    gata3_vep_record = {
        **VEP_SAMPLE_RECORD,
        "input": "NC_000010.11:g.8045419A>G",
        "most_severe_consequence": "missense_variant",
        "transcript_consequences": [
            {
                "transcript_id": "NM_001002295.2",
                "hgvsc": "NM_001002295.2:c.1A>G",
                "hgvsp": "NP_001002295.1:p.Met1Val",
                "mane_select": "ENST00000379328.9",
                "mane": ["MANE_Select"],
                "consequence_terms": ["missense_variant"],
                "protein_start": 1,
                "protein_end": 1,
                "amino_acids": "M/V",
                "cadd_phred": 22.0,
                "cadd_raw": 4.1,
                "fathmm-xf_coding_pred": "D",
                "fathmm-xf_coding_score": 0.8,
                "fathmm-xf_coding_rankscore": 0.7,
                "phylop100way_vertebrate": 4.0,
                "gene_symbol": "GATA3",
            }
        ],
        "colocated_variants": [],
    }

    monkeypatch.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        lambda submitted_variant: gata3_record,
    )
    monkeypatch.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: gata3_vep_record,
    )

    response_without_tumor = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_001002295.2:c.1A>G"},
    )
    response_with_tumor = client.get(
        "/summarizeEvidence",
        params={
            "variant": "NM_001002295.2:c.1A>G",
            "tumorType": "Hodgkin Lymphoma",
        },
    )

    assert response_without_tumor.status_code == 200
    functional_without_tumor = response_without_tumor.json()["oncogenicityEvidence"]["functional"]
    assert functional_without_tumor["status"] == "not_available"

    assert response_with_tumor.status_code == 200
    functional_with_tumor = response_with_tumor.json()["oncogenicityEvidence"]["functional"]
    assert functional_with_tumor["status"] == "applied"
    assert functional_with_tumor["evidenceCode"] == "OS2"


def test_predictive_gata3_uses_tumor_type_for_ovs1(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gata3_record = {
        "@id": "http://reg.genome.network/allele/CATEST3",
        "genomicAlleles": [
            {
                "chromosome": "10",
                "coordinates": [
                    {
                        "allele": "A",
                        "end": 8055658,
                        "referenceAllele": "G",
                        "start": 8055657,
                    }
                ],
                "hgvs": ["NC_000010.11:g.8055658G>A"],
                "referenceGenome": "GRCh38",
            }
        ],
        "transcriptAlleles": [
            {
                "geneSymbol": "GATA3",
                "geneNCBI_id": 2625,
                "hgvs": ["NM_001002295.2:c.3G>A"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "10",
                        "start": 8116731,
                        "end": 8116732,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_001002295.1:p.Met1Ile",
                    "hgvsWellDefined": "NP_001002295.1:p.Met1Ile",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_001002295.2:c.3G>A",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_001002295.1:p.Met1Ile",
                        }
                    },
                },
            }
        ],
    }

    gata3_vep_record = {
        **VEP_SAMPLE_RECORD,
        "input": "NC_000010.11:g.8055658G>A",
        "most_severe_consequence": "start_lost",
        "transcript_consequences": [
            {
                "transcript_id": "NM_001002295.2",
                "hgvsc": "NM_001002295.2:c.3G>A",
                "hgvsp": "NP_001002295.1:p.Met1Ile",
                "mane_select": "ENST00000379328.9",
                "mane": ["MANE_Select"],
                "consequence_terms": ["start_lost"],
                "protein_start": 1,
                "protein_end": 1,
                "amino_acids": "M/I",
                "cadd_phred": 27.8,
                "cadd_raw": 4.97372,
                "fathmm-xf_coding_pred": "N",
                "fathmm-xf_coding_score": 0.099101,
                "fathmm-xf_coding_rankscore": 0.19857,
                "phylop100way_vertebrate": 9.496,
                "gene_symbol": "GATA3",
            }
        ],
        "colocated_variants": [],
    }

    monkeypatch.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        lambda submitted_variant: gata3_record,
    )
    monkeypatch.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: gata3_vep_record,
    )

    response_without_tumor = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_001002295.2:c.3G>A"},
    )
    response_tsg = client.get(
        "/summarizeEvidence",
        params={
            "variant": "NM_001002295.2:c.3G>A",
            "tumorType": "Breast Cancer",
        },
    )
    response_oncogene = client.get(
        "/summarizeEvidence",
        params={
            "variant": "NM_001002295.2:c.3G>A",
            "tumorType": "Hodgkin Lymphoma",
        },
    )

    predictive_without_tumor = response_without_tumor.json()["oncogenicityEvidence"][
        "predictive"
    ]
    assert predictive_without_tumor["status"] == "not_available"
    assert "requires tumor type" in predictive_without_tumor["evidenceStatement"]

    predictive_tsg = response_tsg.json()["oncogenicityEvidence"]["predictive"]
    assert predictive_tsg["status"] == "applied"
    assert predictive_tsg["score"] == 8
    assert predictive_tsg["evidenceCode"] == "OVS1"
    assert predictive_tsg["matchedData"]["resolvedGeneRole"] == "tsg"

    predictive_oncogene = response_oncogene.json()["oncogenicityEvidence"]["predictive"]
    assert predictive_oncogene["status"] == "applied"
    assert predictive_oncogene["score"] == 0
    assert predictive_oncogene["evidenceCode"] is None
    assert predictive_oncogene["matchedData"]["resolvedGeneRole"] == "oncogene"
    assert predictive_oncogene["evidenceStatement"] == (
        "Predictive evidence did not meet current scoring criteria."
    )


def test_functional_conflicting_clinmave_rows_return_applied_zero(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    clinmave_dir = tmp_path / "clinmave"
    clinmave_dir.mkdir()
    (clinmave_dir / "variants.ATM.csv").write_text(
        '"Identifier","Chrom","Position","Ref/Alt","Gene name","Score","Dataset ID","Molecular consequence","Functional description","Phenotype","Cross-assay hits","Functional classification","ClinVar information","Population frequency","TCGA summary","MAVE technique","Mutagenesis strategy","Publication"\n'  # noqa: E501
        '"NM_000051.4(ATM):c.797G>A (p.Trp266Ter)","chr11","108244922","G/A","ATM","-1.0","dataset0001","Nonsense","Reduced ATM activity","ATM phenotype","5","Loss-of-function","","","NA","CRISPR-Based Genome Editing","Base editing","12345678"\n'  # noqa: E501
        '"NM_000051.4(ATM):c.797G>A (p.Trp266Ter)","chr11","108244922","G/A","ATM","0.5","dataset0002","Nonsense","Neutral ATM activity","ATM phenotype","5","Functionally normal","","","NA","CRISPR-Based Genome Editing","Base editing","12345678"\n',  # noqa: E501
        encoding="utf-8",
    )

    monkeypatch.setattr(functional, "CLINMAVE_DATA_DIR", clinmave_dir)
    functional._GENE_INDEX_CACHE.clear()

    atm_record = {
        "@id": "http://reg.genome.network/allele/CATEST2",
        "genomicAlleles": [
            {
                "chromosome": "11",
                "coordinates": [
                    {
                        "allele": "A",
                        "end": 108244922,
                        "referenceAllele": "G",
                        "start": 108244921,
                    }
                ],
                "hgvs": ["NC_000011.10:g.108244922G>A"],
                "referenceGenome": "GRCh38",
            }
        ],
        "transcriptAlleles": [
            {
                "geneSymbol": "ATM",
                "geneNCBI_id": 472,
                "hgvs": ["NM_000051.4:c.797G>A"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "11",
                        "start": 108111262,
                        "end": 108111263,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_000042.3:p.Trp266Ter",
                    "hgvsWellDefined": "NP_000042.3:p.Trp266Ter",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_000051.4:c.797G>A",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_000042.3:p.Trp266Ter",
                        }
                    },
                },
            }
        ],
    }

    atm_vep_record = {
        **VEP_SAMPLE_RECORD,
        "input": "NC_000011.10:g.108244922G>A",
        "most_severe_consequence": "stop_gained",
        "transcript_consequences": [
            {
                "transcript_id": "NM_000051.4",
                "hgvsc": "NM_000051.4:c.797G>A",
                "hgvsp": "NP_000042.3:p.Trp266Ter",
                "mane_select": "ENST00000278616.10",
                "mane": ["MANE_Select"],
                "consequence_terms": ["stop_gained"],
                "protein_start": 266,
                "protein_end": 266,
                "amino_acids": "W/*",
                "cadd_phred": 36.0,
                "cadd_raw": 7.0,
                "fathmm-xf_coding_pred": "D",
                "fathmm-xf_coding_score": 0.9,
                "fathmm-xf_coding_rankscore": 0.9,
                "phylop100way_vertebrate": 5.0,
                "gene_symbol": "ATM",
            }
        ],
        "colocated_variants": [],
    }

    monkeypatch.setattr(
        variant_normalizer,
        "fetch_allele_registry_record",
        lambda submitted_variant: atm_record,
    )
    monkeypatch.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: atm_vep_record,
    )

    response = client.get(
        "/summarizeEvidence",
        params={"variant": "NM_000051.4:c.797G>A"},
    )

    assert response.status_code == 200
    functional_evidence = response.json()["oncogenicityEvidence"]["functional"]
    assert functional_evidence["status"] == "applied"
    assert functional_evidence["score"] == 0
    assert "conflicting functional classifications" in functional_evidence["evidenceStatement"]
    assert functional_evidence["dataAbsentReason"] is None
    assert functional_evidence["matchedData"]["matchingRowCount"] == 2


def test_predict_batch_returns_observations() -> None:
    response = client.post(
        "/predictOncogenicity",
        json={"variants": ["NM_004119.3:c.2073T>G", "ENST00000241453.12:c.2073T>G"]},
    )

    assert response.status_code == 200
    body = response.json()
    assert len(body["observations"]) == 2
    assert body["observations"][0]["valueInteger"] == 2
    assert body["observations"][1]["valueInteger"] == 2


def test_predictive_os1_matches_exact_protein_alias(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        predictive,
        "search_clinvar_variation_ids",
        lambda query, retmax=100: ["40364", "13964"],
    )
    monkeypatch.setattr(
        predictive,
        "fetch_clinvar_summaries",
        lambda variation_ids: {
            "40364": {
                "protein_change": "G464V, G427V, G376V, G442V, G412V, G430V, G467V, G504V",
                "oncogenicity_classification": {"description": "Oncogenic"},
                "title": "NM_004333.6(BRAF):c.1391G>T (p.Gly464Val)",
            },
            "13964": {
                "protein_change": "G464E, G427E, G442E, G504E, G412E, G467E, G376E, G430E",
                "oncogenicity_classification": {"description": "Oncogenic"},
                "title": "NM_004333.6(BRAF):c.1391G>A (p.Gly464Glu)",
            },
        },
    )

    braf_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                "geneSymbol": "BRAF",
                "geneNCBI_id": 673,
                "hgvs": ["NM_004333.6:c.1391G>T"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "7",
                        "start": 140453135,
                        "end": 140453135,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_004324.2:p.Gly464Val",
                    "hgvsWellDefined": "NP_004324.2:p.Gly464Val",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_004333.6:c.1391G>T",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_004324.2:p.Gly464Val",
                        }
                    },
                },
            }
        ],
    }

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: braf_record
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NC_000007.14:g.140781617C>A"},
        )
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert response.status_code == 200
    body = response.json()
    predictive_evidence = body["oncogenicityEvidence"]["predictive"]
    assert predictive_evidence["score"] == 4
    assert predictive_evidence["evidenceCode"] == "OS1"
    assert predictive_evidence["matchedData"]["matchedVariationId"] == "40364"
    assert predictive_evidence["matchedData"]["proteinChangeOneLetter"] == "G464V"


def test_predictive_os1_rejects_same_residue_different_protein_change(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        predictive,
        "search_clinvar_variation_ids",
        lambda query, retmax=100: ["13964"],
    )
    monkeypatch.setattr(
        predictive,
        "fetch_clinvar_summaries",
        lambda variation_ids: {
            "13964": {
                "protein_change": "G464E, G427E, G442E, G504E, G412E, G467E, G376E, G430E",
                "oncogenicity_classification": {"description": "Oncogenic"},
                "title": "NM_004333.6(BRAF):c.1391G>A (p.Gly464Glu)",
            }
        },
    )

    braf_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                "geneSymbol": "BRAF",
                "geneNCBI_id": 673,
                "hgvs": ["NM_004333.6:c.1391G>T"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "7",
                        "start": 140453135,
                        "end": 140453135,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_004324.2:p.Gly464Val",
                    "hgvsWellDefined": "NP_004324.2:p.Gly464Val",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_004333.6:c.1391G>T",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_004324.2:p.Gly464Val",
                        }
                    },
                },
            }
        ],
    }

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: braf_record
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NC_000007.14:g.140781617C>A"},
        )
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert response.status_code == 200
    body = response.json()
    predictive_evidence = body["oncogenicityEvidence"]["predictive"]
    assert predictive_evidence["score"] == 2
    assert predictive_evidence["evidenceCode"] == "OM4"


def test_predictive_om4_matches_different_same_residue_oncogenic_variant(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_search_clinvar_variation_ids(query: str, retmax: int = 100) -> list[str]:
        if query == "BRAF[gene] AND G464V[varnam]":
            return []
        if query == "BRAF[gene] AND Gly464Val[varnam]":
            return []
        if query == "BRAF[gene] AND G464[varnam]":
            return ["13964"]
        if query == "BRAF[gene] AND Gly464[varnam]":
            return ["13964"]
        return []

    monkeypatch.setattr(
        predictive,
        "search_clinvar_variation_ids",
        fake_search_clinvar_variation_ids,
    )
    monkeypatch.setattr(
        predictive,
        "fetch_clinvar_summaries",
        lambda variation_ids: {
            "13964": {
                "protein_change": "G464E, G427E, G442E, G504E, G412E, G467E, G376E, G430E",
                "oncogenicity_classification": {"description": "Likely oncogenic"},
                "title": "NM_004333.6(BRAF):c.1391G>A (p.Gly464Glu)",
            }
        },
    )

    braf_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                "geneSymbol": "BRAF",
                "geneNCBI_id": 673,
                "hgvs": ["NM_004333.6:c.1391G>T"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "7",
                        "start": 140453135,
                        "end": 140453135,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_004324.2:p.Gly464Val",
                    "hgvsWellDefined": "NP_004324.2:p.Gly464Val",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_004333.6:c.1391G>T",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_004324.2:p.Gly464Val",
                        }
                    },
                },
            }
        ],
    }

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: braf_record
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NC_000007.14:g.140781617C>A"},
        )
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert response.status_code == 200
    body = response.json()
    predictive_evidence = body["oncogenicityEvidence"]["predictive"]
    assert predictive_evidence["score"] == 2
    assert predictive_evidence["evidenceCode"] == "OM4"
    assert predictive_evidence["matchedData"]["matchedVariationId"] == "13964"
    assert predictive_evidence["matchedData"]["matchedRule"] == "OM4"


def test_predictive_ovs1_applies_for_frameshift_in_tsg() -> None:
    tsg_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                **CLINGEN_SAMPLE_RECORD["transcriptAlleles"][0],
                "geneSymbol": "TP53",
                "geneNCBI_id": 7157,
                "hgvs": ["NM_000546.6:c.375_376del"],
            }
        ],
    }
    frameshift_record = {
        **VEP_SAMPLE_RECORD,
        "most_severe_consequence": "frameshift_variant",
        "transcript_consequences": [
            {
                **VEP_SAMPLE_RECORD["transcript_consequences"][0],
                "transcript_id": "NM_000546.6",
                "consequence_terms": ["frameshift_variant"],
                "gene_symbol": "TP53",
                "hgvsc": "NM_000546.6:c.375_376del",
                "hgvsp": "NP_000537.3:p.Lys125fs",
            }
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: frameshift_record,
    )
    try:
        original_fetch = variant_normalizer.fetch_allele_registry_record
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: tsg_record
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NM_000546.6:c.375_376del"},
        )
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    predictive_evidence = body["oncogenicityEvidence"]["predictive"]
    assert predictive_evidence["score"] == 8
    assert predictive_evidence["evidenceCode"] == "OVS1"
    assert predictive_evidence["matchedData"]["geneRole"] == "tsg"


def test_predictive_om2_applies_for_inframe_deletion_in_oncogene() -> None:
    inframe_deletion_record = {
        **VEP_SAMPLE_RECORD,
        "most_severe_consequence": "inframe_deletion",
        "transcript_consequences": [
            {
                **VEP_SAMPLE_RECORD["transcript_consequences"][0],
                "consequence_terms": ["inframe_deletion"],
                "protein_start": 691,
                "protein_end": 691,
                "amino_acids": "F/-",
                "hgvsp": "NP_004110.2:p.Phe691del",
            }
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: inframe_deletion_record,
    )
    try:
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NM_004119.3:c.2073_2075del"},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    predictive_evidence = body["oncogenicityEvidence"]["predictive"]
    assert predictive_evidence["score"] == 2
    assert predictive_evidence["evidenceCode"] == "OM2"
    assert predictive_evidence["matchedData"]["geneRole"] == "oncogene"


def test_predictive_sbp2_applies_for_low_phylo_p_synonymous_variant() -> None:
    synonymous_record = {
        **VEP_SAMPLE_RECORD,
        "most_severe_consequence": "synonymous_variant",
        "transcript_consequences": [
            {
                **VEP_SAMPLE_RECORD["transcript_consequences"][0],
                "consequence_terms": ["synonymous_variant"],
                "hgvsp": "NP_004110.2:p.Phe691=",
                "phylop100way_vertebrate": 1.2,
                "cadd_phred": 10.1,
                "cadd_raw": 0.1,
            }
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: synonymous_record,
    )
    try:
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NM_004119.3:c.2073T>C"},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    predictive_evidence = body["oncogenicityEvidence"]["predictive"]
    assert predictive_evidence["score"] == -1
    assert predictive_evidence["evidenceCode"] == "SBP2"
    assert predictive_evidence["matchedData"]["phyloP100wayVertebrate"] == 1.2


def test_predict_single_returns_hotspot_component_when_hotspot_matches() -> None:
    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(hotspots, "_get_hotspot_index", lambda: FAKE_HOTSPOT_INDEX)
    try:
        response = client.get(
            "/predictOncogenicity",
            params={"variant": "NM_004119.3:c.2073T>G"},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    assert body["valueInteger"] == 6

    hotspot_component = _component_by_code(body, "hotspots-evidence")
    assert hotspot_component["valueInteger"] == 4
    assert hotspot_component["interpretation"][0]["coding"][0]["code"] == "OS3"
    assert hotspot_component["interpretation"][0]["text"] == (
        "Located in Cancer Hotspots with at least 50 observed mutations (75) "
        "and at least 10 occurrences of the same protein event (12)."
    )


def test_summarize_evidence_returns_hotspot_match_data_when_present() -> None:
    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(hotspots, "_get_hotspot_index", lambda: FAKE_HOTSPOT_INDEX)
    try:
        response = client.get(
            "/summarizeEvidence",
            params={"variant": "NM_004119.3:c.2073T>G"},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    assert body["overallScore"] == 6
    assert body["oncogenicityEvidence"]["hotspots"]["evidenceCode"] == "OS3"
    assert body["oncogenicityEvidence"]["hotspots"]["matchedData"]["gene"] == "FLT3"
    assert body["oncogenicityEvidence"]["hotspots"]["matchedData"]["position"] == 691
    assert body["oncogenicityEvidence"]["hotspots"]["matchedData"]["proteinHgvs"] == "p.F691L"


def test_hotspot_indel_match_uses_protein_hgvs_bounds_before_vep_positions() -> None:
    indel_index = hotspots.HotspotIndex(
        snv_by_gene={},
        indel_by_gene={
            "EGFR": (
                hotspots.IndelHotspotRecord(
                    gene="EGFR",
                    hotspot_position="745-759",
                    event="E746_A750del",
                    event_start=746,
                    event_end=750,
                    mutation_count=156,
                    variant_count=123,
                ),
            )
        },
    )
    annotated_variant = AnnotatedVariant(
        normalizedVariant=NormalizedVariant(
            submitted_variant="NM_005228.5:c.2235_2249del",
            normalization=NormalizationMetadata(
                source="test",
                queried_variant="NM_005228.5:c.2235_2249del",
            ),
            geneSymbol="EGFR",
        ),
        annotationStatus="complete",
        basicAnnotation=BasicAnnotation(
            mostSevereConsequence="inframe_deletion",
            transcriptConsequences=[
                TranscriptConsequence(
                    transcriptRefSeq="NM_005228.5",
                    consequenceTerms=["inframe_deletion"],
                    proteinStart=745,
                    proteinEnd=750,
                    proteinHgvs="p.E746_A750del",
                    proteinEventType="deletion",
                    isManeSelect=True,
                )
            ],
        ),
        computationalAnnotation=None,
    )

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(hotspots, "_get_hotspot_index", lambda: indel_index)
    try:
        evidence = hotspots.build_hotspots_evidence(annotated_variant)
    finally:
        monkeypatch_context.undo()

    assert evidence.evidenceCode == "OS3"
    assert evidence.score == 4
    assert evidence.matchedData["event"] == "E746_A750del"
    assert evidence.matchedData["proteinStart"] == 745
    assert evidence.matchedData["proteinEnd"] == 750


def test_hotspot_indel_match_still_rejects_event_mismatch() -> None:
    indel_index = hotspots.HotspotIndex(
        snv_by_gene={},
        indel_by_gene={
            "EGFR": (
                hotspots.IndelHotspotRecord(
                    gene="EGFR",
                    hotspot_position="745-759",
                    event="L747_T751del",
                    event_start=747,
                    event_end=751,
                    mutation_count=156,
                    variant_count=9,
                ),
            )
        },
    )
    annotated_variant = AnnotatedVariant(
        normalizedVariant=NormalizedVariant(
            submitted_variant="NM_005228.5:c.2235_2249del",
            normalization=NormalizationMetadata(
                source="test",
                queried_variant="NM_005228.5:c.2235_2249del",
            ),
            geneSymbol="EGFR",
        ),
        annotationStatus="complete",
        basicAnnotation=BasicAnnotation(
            mostSevereConsequence="inframe_deletion",
            transcriptConsequences=[
                TranscriptConsequence(
                    transcriptRefSeq="NM_005228.5",
                    consequenceTerms=["inframe_deletion"],
                    proteinStart=745,
                    proteinEnd=750,
                    proteinHgvs="p.E746_A750del",
                    proteinEventType="deletion",
                    isManeSelect=True,
                )
            ],
        ),
        computationalAnnotation=None,
    )

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(hotspots, "_get_hotspot_index", lambda: indel_index)
    try:
        evidence = hotspots.build_hotspots_evidence(annotated_variant)
    finally:
        monkeypatch_context.undo()

    assert evidence.evidenceCode is None
    assert evidence.score == 0
    assert evidence.evidenceStatement == "Hotspots evidence did not meet current scoring criteria."


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


def test_normalization_converts_deletion_protein_notation_to_one_letter() -> None:
    deletion_record = {
        **CLINGEN_SAMPLE_RECORD,
        "transcriptAlleles": [
            {
                "geneSymbol": "CFTR",
                "geneNCBI_id": 1080,
                "hgvs": ["NM_000492.4:c.1521_1523delCTT"],
                "genomeAlignments": [
                    {
                        "referenceGenome": "GRCh37",
                        "chromosome": "7",
                        "start": 117559592,
                        "end": 117559594,
                    }
                ],
                "proteinEffect": {
                    "hgvs": "NP_000483.3:p.Phe508del",
                    "hgvsWellDefined": "NP_000483.3:p.Phe508del",
                },
                "MANE": {
                    "maneStatus": "MANE Select",
                    "nucleotide": {
                        "RefSeq": {
                            "hgvs": "NM_000492.4:c.1521_1523delCTT",
                        }
                    },
                    "protein": {
                        "RefSeq": {
                            "hgvs": "NP_000483.3:p.Phe508del",
                        }
                    },
                },
            }
        ],
    }

    original_fetch = variant_normalizer.fetch_allele_registry_record
    try:
        variant_normalizer.fetch_allele_registry_record = lambda submitted_variant: deletion_record
        normalized_variant = variant_normalizer.normalize_variant("NM_000492.4:c.1521_1523delCTT")
    finally:
        variant_normalizer.fetch_allele_registry_record = original_fetch

    assert normalized_variant.protein.hgvs_protein_full == "NP_000483.3:p.Phe508del"
    assert normalized_variant.protein.hgvs_3letter == "p.Phe508del"
    assert normalized_variant.protein.hgvs_1letter == "p.F508del"
    assert normalized_variant.protein.short_name == "F508del"
    assert normalized_variant.protein.civic_profile_name == "CFTR F508del"


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
    assert body["valueInteger"] == 0
    population = _component_by_code(body, "population-evidence")
    assert "valueInteger" not in population
    assert population["interpretation"] == []
    assert population["dataAbsentReason"] == {
        "coding": [
            {
                "system": "http://terminology.hl7.org/CodeSystem/data-absent-reason",
                "code": "error",
                "display": "error",
            }
        ],
        "text": "Population evidence could not be evaluated because annotation data was unavailable.",
    }
    computational = _component_by_code(body, "computational-evidence")
    assert computational["dataAbsentReason"] == {
        "coding": [
            {
                "system": "http://terminology.hl7.org/CodeSystem/data-absent-reason",
                "code": "error",
                "display": "error",
            }
        ],
        "text": "Computational evidence could not be evaluated because annotation data was unavailable.",
    }


def test_annotate_single_returns_failed_annotation_payload_when_vep_fails() -> None:
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
            "/annotateVariant",
            params={"variant": "NM_004119.3:c.2073T>G"},
        )
    finally:
        monkeypatch_context.undo()

    assert response.status_code == 200
    body = response.json()
    assert body["normalizedVariant"]["submitted_variant"] == "NM_004119.3:c.2073T>G"
    assert body["annotationStatus"] == "failed"
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
    first_population = _component_by_code(body["observations"][0], "population-evidence")
    second_population = _component_by_code(body["observations"][1], "population-evidence")
    assert first_population["valueInteger"] == 1
    assert second_population["dataAbsentReason"]["coding"][0]["code"] == "error"


def test_predict_single_returns_op4_when_population_data_is_missing() -> None:
    no_population_record = {
        **VEP_SAMPLE_RECORD,
        "colocated_variants": [],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: no_population_record,
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
    assert body["valueInteger"] == 2
    population = _component_by_code(body, "population-evidence")
    assert population["valueInteger"] == 1
    assert population["interpretation"][0]["coding"][0]["code"] == "OP4"
    assert population["interpretation"][0]["text"] == "Absent from gnomAD controls."


def test_predict_single_returns_sbp1_for_low_cadd_and_benign_fathmm_xf() -> None:
    low_cadd_benign_record = {
        **VEP_SAMPLE_RECORD,
        "transcript_consequences": [
            {
                **VEP_SAMPLE_RECORD["transcript_consequences"][0],
                "cadd_phred": 10.4,
                "cadd_raw": 0.42,
                "fathmm-xf_coding_pred": "N",
                "fathmm-xf_coding_score": 0.12,
                "fathmm-xf_coding_rankscore": 0.08,
            },
            VEP_SAMPLE_RECORD["transcript_consequences"][1],
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: low_cadd_benign_record,
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
    assert body["valueInteger"] == 0
    computational = _component_by_code(body, "computational-evidence")
    assert computational["valueInteger"] == -1
    assert computational["interpretation"][0]["coding"][0]["code"] == "SBP1"


def test_predict_single_applies_op1_for_high_cadd_outside_missense() -> None:
    non_missense_record = {
        **VEP_SAMPLE_RECORD,
        "most_severe_consequence": "synonymous_variant",
        "transcript_consequences": [
            {
                **VEP_SAMPLE_RECORD["transcript_consequences"][0],
                "consequence_terms": ["synonymous_variant"],
                "cadd_phred": 25.3,
                "fathmm-xf_coding_pred": "N",
            },
            VEP_SAMPLE_RECORD["transcript_consequences"][1],
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: non_missense_record,
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
    assert body["valueInteger"] == 2
    computational = _component_by_code(body, "computational-evidence")
    assert computational["valueInteger"] == 1
    assert computational["interpretation"][0]["coding"][0]["code"] == "OP1"
    assert computational["interpretation"][0]["text"] == (
        "CADD supports oncogenicity for this variant (PHRED 25.3; most severe consequence synonymous_variant)."
    )


def test_predict_single_returns_no_computational_code_for_low_cadd_non_missense() -> None:
    low_cadd_non_missense_record = {
        **VEP_SAMPLE_RECORD,
        "most_severe_consequence": "inframe_deletion",
        "transcript_consequences": [
            {
                **VEP_SAMPLE_RECORD["transcript_consequences"][0],
                "consequence_terms": ["inframe_deletion"],
                "cadd_phred": 10.4,
                "cadd_raw": 0.42,
                "fathmm-xf_coding_pred": "N",
            },
            VEP_SAMPLE_RECORD["transcript_consequences"][1],
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: low_cadd_non_missense_record,
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
    assert body["valueInteger"] == 3
    computational = _component_by_code(body, "computational-evidence")
    assert computational["valueInteger"] == 0
    assert computational["interpretation"][0]["coding"] == []
    assert computational["interpretation"][0]["text"] == (
        "Computational missense benign rules were not applicable because the most severe consequence "
        "was inframe_deletion."
    )
    predictive = _component_by_code(body, "predictive-evidence")
    assert predictive["valueInteger"] == 2
    assert predictive["interpretation"][0]["coding"][0]["code"] == "OM2"


def test_predict_single_returns_sbvs1_for_high_population_frequency() -> None:
    high_population_record = {
        **VEP_SAMPLE_RECORD,
        "colocated_variants": [
            {
                "frequencies": {
                    "C": {
                        "gnomade_nfe": 0.054,
                        "gnomade": 0.041,
                        "gnomadg": 0.039,
                    }
                }
            }
        ],
    }

    monkeypatch_context = pytest.MonkeyPatch()
    monkeypatch_context.setattr(
        variant_annotator,
        "fetch_vep_annotation_record",
        lambda normalized_variant: high_population_record,
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
    assert body["valueInteger"] == -7
    population = _component_by_code(body, "population-evidence")
    assert population["valueInteger"] == -8
    assert population["interpretation"][0]["coding"][0]["code"] == "SBVS1"


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
