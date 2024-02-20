#!/usr/bin/env python3

import pytest

from variant_classification.acmg_rules.ps3 import Ps3
from variant_classification.acmg_rules.bs3 import Bs3
from variant_classification.acmg_rules.functional_splicing_assay_utils import (
    summarise_func_data,
    adjust_strength_according_to_rna_data_pvs1,
    assess_splicing_data_bp7,
)
from variant_classification.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
)
from variant_classification.variant import FunctionalData, RNAData, ALLELIC

func_no_result = FunctionalData(pathogenic=False, benign=False)
func_both_true = FunctionalData(pathogenic=True, benign=True)
func_benign = FunctionalData(pathogenic=False, benign=True)
func_patho = FunctionalData(pathogenic=True, benign=False)


def test_both_true_error():
    """
    Test that error is raised when both benign and pathogenic are true in FunctionalData
    """
    with pytest.raises(AttributeError):
        out = summarise_func_data([func_both_true])


def test_conflicting_result():
    """
    Test output with confliciting results
    """
    confliciting_results = [func_benign, func_patho]
    result_bs3 = Bs3.assess_rule(confliciting_results)
    result_ps3 = Ps3.assess_rule(confliciting_results)
    assert not result_bs3.status and not result_ps3.status


def test_ben_result():
    """
    Test output for benign result
    """
    benign_results = [func_benign, func_benign]
    result_bs3 = Bs3.assess_rule(benign_results)
    result_ps3 = Ps3.assess_rule(benign_results)
    assert not result_ps3.status and result_bs3.status


def test_patho_result():
    """
    Test output for pathogenic result
    """
    patho_results = [func_patho, func_patho]
    result_bs3 = Bs3.assess_rule(patho_results)
    result_ps3 = Ps3.assess_rule(patho_results)
    assert result_ps3.status and not result_bs3.status


minigene_no_quant = RNAData(True, False, ALLELIC.CONSTRUCT, None)
minigene_quant_patho_strong = RNAData(True, False, ALLELIC.CONSTRUCT, 0.95)
minigene_quant_patho = RNAData(True, False, ALLELIC.CONSTRUCT, 0.85)
minigene_quant_indiff = RNAData(True, False, ALLELIC.CONSTRUCT, 0.75)
minigene_quant_ben = RNAData(True, False, ALLELIC.CONSTRUCT, 0.65)
patrna_no_quant = RNAData(False, True, ALLELIC.TRUE, None)
patrna_no_allelic_ben = RNAData(False, True, ALLELIC.FALSE, 0)
patrna_no_allelic_patho = RNAData(False, True, ALLELIC.FALSE, 0.98)
patrna_no_allelic_indiff = RNAData(False, True, ALLELIC.FALSE, 0.5)
patrna_quant_patho_strong = RNAData(False, True, ALLELIC.TRUE, 0.95)
patrna_quant_patho = RNAData(False, True, ALLELIC.TRUE, 0.85)
patrna_quant_indiff = RNAData(False, True, ALLELIC.TRUE, 0.75)
patrna_quant_benign = RNAData(False, True, ALLELIC.TRUE, 0.65)


def test_empty_df():
    """
    Test empty data frame
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1([], example_result)
    assert (
        adjusted_rule.status
        and adjusted_rule.strength is evidence_strength.VERY_STRONG
        and adjusted_rule.name == "PVS1"
    )


def test_minigene_no_quant():
    """
    Test rule strenght modification for minigene assay without quantification
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [minigene_no_quant], example_result
    )
    assert (
        adjusted_rule.status
        and adjusted_rule.strength is evidence_strength.VERY_STRONG
        and adjusted_rule.name == "PVS1"
    )


def test_mingene_quant_patho_strong_pvs1():
    """
    Test rule strenght modification for minigene with quantification meeting threshold of 0.9
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [minigene_quant_patho_strong], example_result
    )
    assert (
        adjusted_rule.status
        and adjusted_rule.strength.value is evidence_strength.STRONG.value
        and adjusted_rule.name == "PVS1_RNA"
    )


def test_mingene_quant_patho_pvs1():
    """
    Test rule strenght modification for minigene with quantification meeting threshold betwenn 0.9 and 0.8
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [minigene_quant_patho], example_result
    )
    assert (
        adjusted_rule.status
        and adjusted_rule.strength.value is evidence_strength.MODERATE.value
        and adjusted_rule.name == "PVS1_RNA"
    )


def test_mingene_quant_ben_pvs1():
    """
    Test rule strenght modification for minigene with quantification meeting threshold between 0.9 and 0.8
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [minigene_quant_ben], example_result
    )
    assert not adjusted_rule.status and adjusted_rule.name == "PVS1_RNA"


def test_mingene_quant_inconclusive_minigene_pvs1():
    """
    Test rule strenght modification for minigene with quantification meeting threshold between 0.9 and 0.8
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [minigene_quant_indiff], example_result
    )
    assert not adjusted_rule.status and adjusted_rule.name == "PVS1_RNA"


def test_patrna_quant_inconclusive_pvs1():
    """
    Test rule strenght modification for minigene with quantification meeting threshold between 0.9 and 0.8
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [patrna_quant_indiff], example_result
    )
    assert not adjusted_rule.status and adjusted_rule.name == "PVS1_RNA"


def test_patrna_benign_pvs1():
    """
    Test rule strenght modification for minigene with quantification meeting threshold between 0.9 and 0.8
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [patrna_quant_benign], example_result
    )
    assert not adjusted_rule.status and adjusted_rule.name == "PVS1_RNA"


def test_pat_rna_quant_patho_strong_pvs1():
    """
    Test rule strength modification for allele specific patient RNA assay with quantification meeting threshold above 0.9
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [patrna_quant_patho_strong], example_result
    )
    assert (
        adjusted_rule.status
        and adjusted_rule.name == "PVS1_RNA"
        and adjusted_rule.strength.value == evidence_strength.VERY_STRONG.value
    )


def test_pat_rna_quant_patho_pvs1():
    """
    Test rule strength modification for allele specific patient RNA assay with quantification meeting threshold between 0.9 and 0.8
    """
    example_result = RuleResult(
        "PVS1",
        rule_type.SPLICING,
        evidence_type.PATHOGENIC,
        True,
        evidence_strength.VERY_STRONG,
        "Example comment.",
    )
    adjusted_rule = adjust_strength_according_to_rna_data_pvs1(
        [patrna_quant_patho], example_result
    )
    assert (
        adjusted_rule.status
        and adjusted_rule.name == "PVS1_RNA"
        and adjusted_rule.strength.value == evidence_strength.STRONG.value
    )


def test_bp7_no_quantification_minigene():
    """
    Test BP7 for no quantification in minigene assay
    """
    performed, result, comment = assess_splicing_data_bp7([minigene_no_quant])
    assert performed and not result


def test_bp7_no_quantification_patrna():
    """
    Test BP7 for no quantification in patient RNA
    """
    performed, result, comment = assess_splicing_data_bp7([patrna_no_quant])
    assert performed and not result


def test_bp7_conflicting_results():
    """
    Test BP7 for conflicting results
    """
    performed, result, comment = assess_splicing_data_bp7(
        [patrna_quant_patho_strong, minigene_quant_ben]
    )
    assert performed and not result


def test_bp7_mingene_benign():
    """
    Test BP7 for minigene benign result
    """
    performed, result, comment = assess_splicing_data_bp7([minigene_quant_ben])
    assert performed and result


def test_bp7_pat_rna_benign():
    """
    Test BP7 for patient RNA benign result
    """
    performed, result, comment = assess_splicing_data_bp7([patrna_quant_benign])
    assert performed and result


def test_bp7_pat_rna_benign_no_allelic():
    """
    Test BP7 for patient RNA benign result
    """
    performed, result, comment = assess_splicing_data_bp7([patrna_no_allelic_ben])
    assert performed and result


def test_bp7_pat_rna_not_benign_no_allelic():
    """
    Test BP7 for patient RNA benign result
    """
    performed, result, comment = assess_splicing_data_bp7([patrna_no_allelic_indiff])
    assert performed and not result
