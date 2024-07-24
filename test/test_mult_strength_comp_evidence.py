#!/usr/bin/env python3

from variant_classification.classify import (
    load_config,
    load_variant,
    get_gene_specific_config,
    check_disease_relevant_transcript,
)
from variant_classification.acmg_rules.computation_evidence_utils import (
    THRESHOLD_DIRECTION,
    Threshold,
)
from variant_classification.acmg_rules.bp4_mult_strength import (
    Bp4_protein_mult_strength,
    Bp4_splicing_mult_strength,
)
from variant_classification.acmg_rules.pp3_mult_strength import (
    Pp3_protein_mult_strength,
    Pp3_splicing_mult_strength,
)
from variant_classification.acmg_rules.utils import evidence_strength
from variant_classification.variant import Variant

import test.paths as paths
from test.test_import_variant import create_json_string_from_variant


splice_thr_patho = Threshold(
    name="SpliceAI",
    direction=THRESHOLD_DIRECTION.GREATER_THAN_OR_EQUAL,
    thresholds=[0.2],
    strengths=[evidence_strength.SUPPORTING],
)

splice_thr_ben = Threshold(
    name="SpliceAI",
    direction=THRESHOLD_DIRECTION.LESS_THAN_OR_EQUAL,
    thresholds=[0.1],
    strengths=[evidence_strength.SUPPORTING],
)

patho_thr_ben = Threshold(
    name="REVEL",
    direction=THRESHOLD_DIRECTION.LESS_THAN_OR_EQUAL,
    thresholds=[0.003, 0.016, 0.183, 0.29],
    strengths=[
        evidence_strength.SUPPORTING,
        evidence_strength.MODERATE,
        evidence_strength.STRONG,
        evidence_strength.VERY_STRONG,
    ],
)

patho_thr_patho = Threshold(
    name="REVEL",
    direction=THRESHOLD_DIRECTION.GREATER_THAN_OR_EQUAL,
    thresholds=[0.932, 0.773, 0.644],
    strengths=[
        evidence_strength.SUPPORTING,
        evidence_strength.MODERATE,
        evidence_strength.STRONG,
    ],
)


def test_splice_supp_ben():
    """
    Test splicing for benign result
    """
    splice_ben = {"SpiceAI": 0.05}
    pp3_result = Pp3_splicing_mult_strength.assess_rule(splice_ben, splice_thr_patho)
    bp4_result = Bp4_splicing_mult_strength.assess_rule(splice_ben, splice_thr_ben)
    assert not pp3_result.status and bp4_result


def test_splice_supp_inconc():
    """
    Test splicing for inconcusive result
    """
    splice_un = {"SpiceAI": 0.15}
    pp3_result = Pp3_splicing_mult_strength.assess_rule(splice_un, splice_thr_patho)
    bp4_result = Bp4_splicing_mult_strength.assess_rule(splice_un, splice_thr_ben)
    assert not pp3_result.status and bp4_result


def test_splice_supp_patho():
    """
    Test splicing for inconcusive result
    """
    splice_patho = {"SpiceAI": 0.25}
    pp3_result = Pp3_splicing_mult_strength.assess_rule(splice_patho, splice_thr_patho)
    bp4_result = Bp4_splicing_mult_strength.assess_rule(splice_patho, splice_thr_ben)
    assert not pp3_result.status and bp4_result


def test_patho_sup_patho():
    """
    Test pathogenicity for supporting pathogenic evidence
    """
    patho_patho_supporting = {"REVEL": 0.7}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info,
        variant.variant_info,
        patho_patho_supporting,
        patho_thr_patho,
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(
        patho_patho_supporting, patho_thr_ben
    )
    assert (
        not bp4_result.status
        and pp3_result.status
        and pp3_result.strength.value == evidence_strength.SUPPORTING.value
    )


def test_patho_mod_patho():
    """
    Test pathogenicity for moderate pathogenic evidence
    """
    patho_patho_moderate = {"REVEL": 0.8}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info,
        variant.variant_info,
        patho_patho_moderate,
        patho_thr_patho,
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(
        patho_patho_moderate, patho_thr_ben
    )
    assert (
        not bp4_result.status
        and pp3_result.status
        and pp3_result.strength.value == evidence_strength.MODERATE.value
    )


def test_patho_strong_patho():
    """
    Test pathogenicity for strong pathogenic evidence
    """
    patho_patho_strong = {"REVEL": 0.98}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info,
        variant.variant_info,
        patho_patho_strong,
        patho_thr_patho,
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(
        patho_patho_strong, patho_thr_ben
    )
    assert (
        not bp4_result.status
        and pp3_result.status
        and pp3_result.strength.value == evidence_strength.STRONG.value
    )


def test_patho_sup_ben():
    """
    Test pathogenicity for supporting benign evidence
    """
    patho_ben_supporting = {"REVEL": 0.2}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info,
        variant.variant_info,
        patho_ben_supporting,
        patho_thr_patho,
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(
        patho_ben_supporting, patho_thr_ben
    )
    assert (
        not pp3_result.status
        and bp4_result.status
        and bp4_result.strength.value == evidence_strength.SUPPORTING.value
    )


def test_patho_mod_ben():
    """
    Test pathogenicity for moderate benign evidence
    """
    patho_ben_moderate = {"REVEL": 0.1}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info,
        variant.variant_info,
        patho_ben_moderate,
        patho_thr_patho,
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(
        patho_ben_moderate, patho_thr_ben
    )
    assert (
        not pp3_result.status
        and bp4_result.status
        and bp4_result.strength.value == evidence_strength.MODERATE.value
    )


def test_patho_strong_ben():
    """
    Test pathogenicity for strong benign evidence
    """
    patho_ben_strong = {"REVEL": 0.01}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info, variant.variant_info, patho_ben_strong, patho_thr_patho
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(patho_ben_strong, patho_thr_ben)
    assert (
        not pp3_result.status
        and bp4_result.status
        and bp4_result.strength.value == evidence_strength.STRONG.value
    )


def test_patho_very_strong_ben():
    """
    Test pathogenicity for very strong benign evidence
    """
    patho_ben_very_strong = {"REVEL": 0.001}
    variant = get_disease_relevant_variant("missense")
    pp3_result = Pp3_protein_mult_strength.assess_rule(
        variant.transcript_info,
        variant.variant_info,
        patho_ben_very_strong,
        patho_thr_patho,
    )
    bp4_result = Bp4_protein_mult_strength.assess_rule(
        patho_ben_very_strong, patho_thr_ben
    )
    assert (
        not pp3_result.status
        and bp4_result.status
        and bp4_result.strength.value == evidence_strength.VERY_STRONG.value
    )


def get_disease_relevant_variant(type: str) -> Variant:
    """
    Create variant object
    """
    # Check which variant is needed for classification
    if type == "splicing":
        path_variant = (
            paths.TEST
            / "test_variants_gene_specific"
            / "ATM_splice_acceptor_variant.json"
        )
    elif type == "missense":
        path_variant = (
            paths.TEST / "test_variants_gene_specific" / "BRCA1_missense_variant.json"
        )
    else:
        raise ValueError(
            f"Type should be either 'splicing' or 'missense' but type is {type}."
        )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_atm.yaml"
    variant = load_variant(var_str)
    config = load_config(path_config)
    final_config = get_gene_specific_config(config, variant.variant_info.gene_name)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    return variant_disease_relevant
