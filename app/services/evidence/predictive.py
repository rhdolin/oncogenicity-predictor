"""Predictive evidence scoring for oncogenicity rules.

The predictive pipeline combines local gene-role heuristics with ClinVar-backed
protein-change matching. ClinVar searches are fielded by gene and variant name,
but because Entrez fielded searches are not exact string lookups, returned
summaries are filtered locally against exact protein aliases.
"""

import os
import re

import httpx

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


CLINVAR_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_TOOL = "oncogenicity-predictor"
NULL_VARIANT_TERMS = {
    "transcript_ablation",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_gained",
}
LENGTH_CHANGE_TERMS = {
    "inframe_insertion",
    "inframe_deletion",
}

OVS1_SCORE = 8
OS1_SCORE = 4
OM2_SCORE = 2
OM4_SCORE = 2
SBP2_SCORE = -1
SBP2_PHYLOP_THRESHOLD = 2.0


class ClinVarLookupError(Exception):
    """Raised when ClinVar cannot be queried or returns an unusable payload."""


def build_predictive_evidence(
    annotated_variant: AnnotatedVariant,
    tumor_type: TumorType | None = None,
) -> EvidenceResult:
    """Evaluate the currently implemented predictive evidence rule set."""
    if (
        annotated_variant.annotationStatus != "complete"
        or annotated_variant.basicAnnotation is None
        or annotated_variant.computationalAnnotation is None
    ):
        return _not_available_result(
            "Predictive evidence could not be evaluated because annotation data was unavailable.",
            matched_data=None,
            reason="error",
        )

    normalized_variant = annotated_variant.normalizedVariant
    consequence = annotated_variant.basicAnnotation.mostSevereConsequence
    phylo_p = annotated_variant.computationalAnnotation.phyloP100wayVertebrate
    gene = normalized_variant.geneSymbol
    gene_role = classify_base_gene_role(gene)
    resolved_gene_role = resolve_gene_role(gene, tumor_type)
    one_letter_change = normalized_variant.protein.short_name
    three_letter_change = _extract_three_letter_protein_token(
        normalized_variant.protein.hgvs_3letter
    )
    one_letter_residue = _extract_one_letter_residue_token(one_letter_change)
    three_letter_residue = _extract_three_letter_residue_token(
        normalized_variant.protein.hgvs_3letter
    )

    matched_data = {
        "implementedRules": ["OVS1", "OS1", "OM2", "SBP2", "OM4"],
        "gene": gene,
        "geneRole": gene_role,
        "resolvedGeneRole": resolved_gene_role,
        "tumorType": tumor_type,
        "mostSevereConsequence": consequence,
        "phyloP100wayVertebrate": phylo_p,
        "proteinChangeOneLetter": one_letter_change,
        "proteinChangeThreeLetter": three_letter_change,
        "residueTokenOneLetter": one_letter_residue,
        "residueTokenThreeLetter": three_letter_residue,
        "attemptedQueries": [],
    }
    blockers: list[tuple[str, str]] = []

    if consequence in NULL_VARIANT_TERMS:
        if not gene:
            blockers.append(
                (
                    (
                        "Predictive OVS1 evidence could not be evaluated because the normalized "
                        "gene symbol was unavailable."
                    ),
                    "unknown",
                )
            )
        elif gene_role == GENE_ROLE_BOTH and resolved_gene_role is None:
            blockers.append(
                (
                    (
                        "Predictive OVS1 evidence could not be evaluated because "
                        f"gene {gene} requires tumor type to resolve predictive interpretation."
                    ),
                    "unsupported",
                )
            )
        elif resolved_gene_role == GENE_ROLE_TSG:
            return EvidenceResult(
                score=OVS1_SCORE,
                evidenceCode="OVS1",
                evidenceStatement=(
                    f"OVS1: {consequence} variant in known tumor suppressor gene {gene}."
                ),
                status="applied",
                source="predictive",
                matchedData=matched_data,
                dataAbsentReason=None,
            )

    if _has_queryable_protein_change(one_letter_change):
        if not gene:
            blockers.append(
                (
                    (
                        "Predictive OS1 evidence could not be evaluated because the normalized "
                        "gene symbol was unavailable."
                    ),
                    "unknown",
                )
            )
        else:
            os1_result = _evaluate_os1(
                gene,
                one_letter_change,
                three_letter_change,
                matched_data,
            )
            if os1_result is not None:
                return os1_result
            if matched_data.get("os1LookupFailed"):
                blockers.append(
                    (
                        "Predictive OS1 evidence could not be evaluated because the ClinVar lookup failed.",
                        "error",
                    )
                )

    if consequence in LENGTH_CHANGE_TERMS or consequence == "stop_lost":
        if not gene:
            blockers.append(
                (
                    (
                        "Predictive OM2 evidence could not be evaluated because the normalized "
                        "gene symbol was unavailable."
                    ),
                    "unknown",
                )
            )
        elif _matches_om2(consequence, gene_role, resolved_gene_role):
            return EvidenceResult(
                score=OM2_SCORE,
                evidenceCode="OM2",
                evidenceStatement=_build_om2_statement(
                    consequence,
                    gene,
                    gene_role,
                    resolved_gene_role,
                ),
                status="applied",
                source="predictive",
                matchedData=matched_data,
                dataAbsentReason=None,
            )

    if consequence == "synonymous_variant":
        if phylo_p is None:
            blockers.append(
                (
                    "Predictive SBP2 evidence could not be evaluated because phyloP annotation was unavailable.",
                    "unknown",
                )
            )
        elif phylo_p < SBP2_PHYLOP_THRESHOLD:
            return EvidenceResult(
                score=SBP2_SCORE,
                evidenceCode="SBP2",
                evidenceStatement=(
                    f"SBP2: Synonymous variant with low conservation scores in {gene or 'unknown gene'}."
                ),
                status="applied",
                source="predictive",
                matchedData=matched_data,
                dataAbsentReason=None,
            )

    if consequence == "missense_variant":
        if not gene or not one_letter_residue:
            blockers.append(
                (
                    (
                        "Predictive OM4 evidence could not be evaluated because the normalized "
                        "residue token was unavailable."
                    ),
                    "unknown",
                )
            )
        else:
            om4_result = _evaluate_om4(
                gene,
                one_letter_change,
                one_letter_residue,
                three_letter_residue,
                matched_data,
            )
            if om4_result is not None:
                return om4_result
            if matched_data.get("om4LookupFailed"):
                blockers.append(
                    (
                        "Predictive OM4 evidence could not be evaluated because the ClinVar lookup failed.",
                        "error",
                    )
                )

    if blockers:
        statement, reason = blockers[0]
        return _not_available_result(statement, matched_data=matched_data, reason=reason)

    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement=(
            "Predictive evidence did not meet current scoring criteria."
        ),
        status="applied",
        source="predictive",
        matchedData=matched_data,
        dataAbsentReason=None,
    )


def _evaluate_os1(
    gene: str,
    one_letter_change: str,
    three_letter_change: str | None,
    matched_data: dict,
) -> EvidenceResult | None:
    queries = _build_os1_queries(gene, one_letter_change, three_letter_change)
    matched_data["os1Queries"] = queries
    matched_data["attemptedQueries"].extend(queries)

    try:
        clinvar_ids, successful_query = _search_first_matching_query(queries)
        matched_data["os1SuccessfulQuery"] = successful_query
        matched_data["os1CandidateVariationIds"] = clinvar_ids
        if not clinvar_ids:
            return None
        summaries = fetch_clinvar_summaries(clinvar_ids)
    except ClinVarLookupError:
        matched_data["os1LookupFailed"] = True
        return None

    match = _find_os1_match(summaries, one_letter_change)
    if match is None:
        return None

    matched_data.update(
        {
            "matchedVariationId": match["uid"],
            "matchedProteinAliases": match["proteinAliases"],
            "matchedOncogenicityClassification": match["oncogenicityClassification"],
            "matchedTitle": match["title"],
            "matchedRule": "OS1",
        }
    )
    return EvidenceResult(
        score=OS1_SCORE,
        evidenceCode="OS1",
        evidenceStatement=(
            "OS1: Same amino acid change as a previously established somatic oncogenic ClinVar variant."
        ),
        status="applied",
        source="predictive",
        matchedData=matched_data,
        dataAbsentReason=None,
    )


def _evaluate_om4(
    gene: str,
    one_letter_change: str | None,
    one_letter_residue: str,
    three_letter_residue: str | None,
    matched_data: dict,
) -> EvidenceResult | None:
    queries = _build_om4_queries(gene, one_letter_residue, three_letter_residue)
    matched_data["om4Queries"] = queries
    matched_data["attemptedQueries"].extend(queries)

    try:
        clinvar_ids, successful_query = _search_first_matching_query(queries)
        matched_data["om4SuccessfulQuery"] = successful_query
        matched_data["om4CandidateVariationIds"] = clinvar_ids
        if not clinvar_ids:
            return None
        summaries = fetch_clinvar_summaries(clinvar_ids)
    except ClinVarLookupError:
        matched_data["om4LookupFailed"] = True
        return None

    match = _find_om4_match(summaries, one_letter_change)
    if match is None:
        return None

    matched_data.update(
        {
            "matchedVariationId": match["uid"],
            "matchedProteinAliases": match["proteinAliases"],
            "matchedOncogenicityClassification": match["oncogenicityClassification"],
            "matchedTitle": match["title"],
            "matchedRule": "OM4",
        }
    )
    return EvidenceResult(
        score=OM4_SCORE,
        evidenceCode="OM4",
        evidenceStatement=(
            "OM4: Missense variant at an amino acid residue where a different "
            "somatic oncogenic missense variant is established in ClinVar."
        ),
        status="applied",
        source="predictive",
        matchedData=matched_data,
        dataAbsentReason=None,
    )


def search_clinvar_variation_ids(query: str, retmax: int = 100) -> list[str]:
    """Return ClinVar variation IDs for a fielded ESearch query."""
    params = {
        "db": "clinvar",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "tool": NCBI_TOOL,
    }
    email = os.getenv("NCBI_EMAIL")
    api_key = os.getenv("NCBI_API_KEY")
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    try:
        response = httpx.get(CLINVAR_ESEARCH_URL, params=params, timeout=30.0)
        response.raise_for_status()
    except httpx.HTTPError as exc:
        raise ClinVarLookupError("ClinVar ESearch failed.") from exc

    payload = response.json()
    results = payload.get("esearchresult", {})
    id_list = results.get("idlist")
    if not isinstance(id_list, list):
        raise ClinVarLookupError("ClinVar ESearch returned an unexpected payload.")
    return [str(uid) for uid in id_list]


def fetch_clinvar_summaries(variation_ids: list[str]) -> dict[str, dict]:
    """Fetch ClinVar ESummary documents keyed by variation ID."""
    if not variation_ids:
        return {}

    params = {
        "db": "clinvar",
        "id": ",".join(variation_ids),
        "retmode": "json",
        "tool": NCBI_TOOL,
    }
    email = os.getenv("NCBI_EMAIL")
    api_key = os.getenv("NCBI_API_KEY")
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    try:
        response = httpx.get(CLINVAR_ESUMMARY_URL, params=params, timeout=30.0)
        response.raise_for_status()
    except httpx.HTTPError as exc:
        raise ClinVarLookupError("ClinVar ESummary failed.") from exc

    payload = response.json()
    result = payload.get("result")
    if not isinstance(result, dict):
        raise ClinVarLookupError("ClinVar ESummary returned an unexpected payload.")
    return {
        str(uid): summary
        for uid, summary in result.items()
        if uid != "uids" and isinstance(summary, dict)
    }


def _build_os1_queries(
    gene: str,
    one_letter_change: str,
    three_letter_change: str | None,
) -> list[str]:
    queries = [f"{gene}[gene] AND {one_letter_change}[varnam]"]
    if three_letter_change and three_letter_change != one_letter_change:
        queries.append(f"{gene}[gene] AND {three_letter_change}[varnam]")
    return queries


def _build_om4_queries(
    gene: str,
    one_letter_residue: str,
    three_letter_residue: str | None,
) -> list[str]:
    queries = [f"{gene}[gene] AND {one_letter_residue}[varnam]"]
    if three_letter_residue and three_letter_residue != one_letter_residue:
        queries.append(f"{gene}[gene] AND {three_letter_residue}[varnam]")
    return queries


def _search_first_matching_query(queries: list[str]) -> tuple[list[str], str | None]:
    for query in queries:
        ids = search_clinvar_variation_ids(query)
        if ids:
            return ids, query
    return [], None


def _extract_three_letter_protein_token(hgvs_3letter: str | None) -> str | None:
    if not hgvs_3letter:
        return None
    return hgvs_3letter.removeprefix("p.") or None


def _extract_one_letter_residue_token(one_letter_change: str | None) -> str | None:
    parsed_change = _parse_one_letter_substitution(one_letter_change)
    if parsed_change is None:
        return None
    return f"{parsed_change['ref']}{parsed_change['position']}"


def _extract_three_letter_residue_token(hgvs_3letter: str | None) -> str | None:
    if not hgvs_3letter:
        return None
    match = re.match(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$", hgvs_3letter)
    if not match:
        return None
    return f"{match.group(1)}{match.group(2)}"


def _find_os1_match(summaries: dict[str, dict], query_change: str) -> dict | None:
    for uid, summary in summaries.items():
        oncogenicity = _extract_oncogenicity_description(summary)
        if not _is_supported_somatic_oncogenicity(oncogenicity):
            continue

        protein_aliases = _extract_protein_aliases(summary.get("protein_change"))
        if query_change not in protein_aliases:
            continue

        return {
            "uid": uid,
            "oncogenicityClassification": oncogenicity,
            "proteinAliases": protein_aliases,
            "title": summary.get("title"),
        }

    return None


def _find_om4_match(summaries: dict[str, dict], query_change: str | None) -> dict | None:
    query_substitution = _parse_one_letter_substitution(query_change)
    if query_substitution is None:
        return None

    for uid, summary in summaries.items():
        oncogenicity = _extract_oncogenicity_description(summary)
        if not _is_supported_somatic_oncogenicity(oncogenicity):
            continue

        protein_aliases = _extract_protein_aliases(summary.get("protein_change"))
        for alias in protein_aliases:
            candidate_substitution = _parse_one_letter_substitution(alias)
            if candidate_substitution is None:
                continue
            if (
                candidate_substitution["ref"] == query_substitution["ref"]
                and candidate_substitution["position"] == query_substitution["position"]
                and candidate_substitution["alt"] != query_substitution["alt"]
            ):
                return {
                    "uid": uid,
                    "oncogenicityClassification": oncogenicity,
                    "proteinAliases": protein_aliases,
                    "title": summary.get("title"),
                }

    return None


def _extract_oncogenicity_description(summary: dict) -> str:
    oncogenicity = summary.get("oncogenicity_classification") or {}
    if not isinstance(oncogenicity, dict):
        return ""
    description = oncogenicity.get("description")
    return description.strip().lower() if isinstance(description, str) else ""


def _is_supported_somatic_oncogenicity(description: str) -> bool:
    return description in {"oncogenic", "likely oncogenic", "oncogenic/likely oncogenic"}


def _extract_protein_aliases(protein_change: object) -> list[str]:
    if not isinstance(protein_change, str):
        return []
    aliases = []
    for token in protein_change.split(","):
        alias = token.strip()
        if alias and alias not in aliases:
            aliases.append(alias)
    return aliases


def _has_queryable_protein_change(one_letter_change: str | None) -> bool:
    if not one_letter_change:
        return False
    return one_letter_change not in {"=", "?"}


def _parse_one_letter_substitution(one_letter_change: str | None) -> dict[str, str] | None:
    if not one_letter_change:
        return None
    match = re.match(r"^([A-Z*])(\d+)([A-Z*])$", one_letter_change)
    if not match:
        return None
    return {
        "ref": match.group(1),
        "position": match.group(2),
        "alt": match.group(3),
    }


def _matches_om2(
    consequence: str | None,
    gene_role: str | None,
    resolved_gene_role: str | None,
) -> bool:
    if consequence in LENGTH_CHANGE_TERMS:
        return gene_role in {GENE_ROLE_ONCOGENE, GENE_ROLE_TSG, GENE_ROLE_BOTH}
    if consequence == "stop_lost":
        if gene_role == GENE_ROLE_BOTH and resolved_gene_role is None:
            return False
        return resolved_gene_role == GENE_ROLE_TSG
    return False


def _build_om2_statement(
    consequence: str | None,
    gene: str,
    gene_role: str | None,
    resolved_gene_role: str | None,
) -> str:
    if consequence == "stop_lost":
        return f"OM2: stop_lost variant in known tumor suppressor gene {gene}."
    if gene_role in {GENE_ROLE_ONCOGENE, GENE_ROLE_BOTH}:
        return f"OM2: {consequence} variant in known cancer gene {gene}."
    if resolved_gene_role == GENE_ROLE_TSG:
        return f"OM2: {consequence} variant in known tumor suppressor gene {gene}."
    return f"OM2: {consequence} variant in known tumor suppressor gene {gene}."


def _not_available_result(
    statement: str,
    matched_data: dict | None,
    reason: str,
) -> EvidenceResult:
    return EvidenceResult(
        score=0,
        evidenceCode=None,
        evidenceStatement=statement,
        status="not_available",
        source="predictive",
        matchedData=matched_data,
        dataAbsentReason=reason,
    )
