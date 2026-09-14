from fastapi.testclient import TestClient
import pytest
from datetime import datetime

from app.main import app
from app.models.annotated_variant import AnnotatedVariant, BasicAnnotation, TranscriptConsequence
from app.models.normalized_variant import NormalizationMetadata, NormalizedVariant
from app.services.annotation import variant_annotator
from app.services.evidence import hotspots
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
    assert predictive["dataAbsentReason"]["coding"][0]["code"] == "unsupported"

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
    assert body["observations"][1]["extension"][0]["valueString"] == "ENST00000241453.12:c.2073T>G"


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
    assert body["valueInteger"] == 1
    computational = _component_by_code(body, "computational-evidence")
    assert computational["valueInteger"] == 0
    assert computational["interpretation"][0]["coding"] == []
    assert computational["interpretation"][0]["text"] == (
        "Computational missense benign rules were not applicable because the most severe consequence "
        "was inframe_deletion."
    )


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
