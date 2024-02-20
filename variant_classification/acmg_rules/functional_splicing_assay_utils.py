#!/usr/bin/env python3

import logging
from typing import Optional

from variant import ALLELIC, FunctionalData, RNAData
from acmg_rules.utils import RuleResult, evidence_strength

logger = logging.getLogger("GenOtoScope_Classify.BS3")


def summarise_func_data(fun_data: list[FunctionalData]) -> tuple[int, int, int]:
    """
    Summarise function data results
    """
    pathogenic_count = 0
    benign_count = 0
    uncertain_count = 0
    for entry in fun_data:
        if entry.pathogenic and entry.benign:
            raise AttributeError(
                f"The result of a functional assay is set to benign and pathogenic. This is not allowed. Please check the input data."
            )
        elif entry.pathogenic:
            pathogenic_count += 1
        elif entry.benign:
            benign_count += 1
        else:
            uncertain_count += 1
    return pathogenic_count, benign_count, uncertain_count


def adjust_strength_according_to_rna_data_pvs1(
    rna_data: list[RNAData], result: RuleResult
) -> RuleResult:
    """
    Modify strength of PVS1 according to available RNA data
    """
    no_quantification = 0
    minigene_90 = False
    minigene = False
    patient_90 = False
    patient = False
    patient_no_allele_specific_patho = False
    patient_no_allele_specific_ben = False
    inconclusive = False
    benign = False
    for entry in rna_data:
        if not entry.quantification or entry.quantification is None:
            no_quantification += 1
        elif entry.minigene and entry.allelic.value is ALLELIC.CONSTRUCT.value:
            if entry.quantification >= 0.9:
                minigene_90 = True
            elif entry.quantification < 0.9 and entry.quantification > 0.8:
                minigene = True
            elif entry.quantification < 0.7:
                benign = True
            else:
                inconclusive = True
        elif entry.patient_rna:
            if entry.allelic.value is ALLELIC.TRUE.value:
                if entry.quantification >= 0.9:
                    patient_90 = True
                elif entry.quantification < 0.9 and entry.quantification > 0.8:
                    patient = True
                elif entry.quantification < 0.7:
                    benign = True
                else:
                    inconclusive = True
            elif entry.allelic.value is ALLELIC.FALSE.value:
                if entry.quantification >= 0.45:
                    patient_no_allele_specific_patho = True
                elif entry.quantification <= 0.05:
                    patient_no_allele_specific_ben = True
                else:
                    inconclusive = True
    strength_adjustment_by_one = {
        evidence_strength.VERY_STRONG.value: evidence_strength.STRONG,
        evidence_strength.STRONG.value: evidence_strength.MODERATE,
        evidence_strength.MODERATE.value: evidence_strength.SUPPORTING,
        evidence_strength.SUPPORTING.value: evidence_strength.SUPPORTING,
    }
    strength_adjustment_by_two = {
        evidence_strength.VERY_STRONG.value: evidence_strength.MODERATE,
        evidence_strength.STRONG.value: evidence_strength.SUPPORTING,
        evidence_strength.MODERATE.value: evidence_strength.SUPPORTING,
        evidence_strength.SUPPORTING.value: evidence_strength.SUPPORTING,
    }
    strength_adjustment_no_allelic_quant = {
        evidence_strength.VERY_STRONG.value: evidence_strength.STRONG,
        evidence_strength.STRONG.value: evidence_strength.MODERATE,
        evidence_strength.MODERATE.value: evidence_strength.MODERATE,
        evidence_strength.SUPPORTING.value: evidence_strength.MODERATE,
    }
    if not rna_data:
        new_comment = "No RNA assay performed."
        update_rule_result(result, comment=new_comment)
        return result
    if no_quantification == len(rna_data):
        new_comment = "No quantification available for the RNA assay. Please check results of RNA assay manually."
        update_rule_result(result, comment=new_comment)
        return result
    result.name = "PVS1_RNA"
    if inconclusive:
        new_comment = "The RNA assay results are inconclusive. Please check results of RNA assay manually."
        update_rule_result(result, status=False, comment=new_comment)
        return result
    if (benign or patient_no_allele_specific_ben) and (
        minigene
        or patient
        or minigene_90
        or patient_90
        or patient_no_allele_specific_patho
    ):
        new_comment = "The RNA assays show contradicting evidence. PVS1_RNA does not apply. Please check resulst of RNA assay manually."
        update_rule_result(result, status=False, comment=new_comment)
        return result
    if benign or patient_no_allele_specific_ben:
        new_comment = "A splicing assay shows no effect on splicing for this variant. PVS1_RNA does therfore not apply."
        update_rule_result(result, status=False, comment=new_comment)
        return result
    if patient_90 and not patient:
        new_comment = "A splicing assay was performed with allele_specific quantification using patient RNA showing >=90% proportion of non-functional transcript."
        update_rule_result(result, comment=new_comment)
        return result
    if patient:
        strength = strength_adjustment_by_one[result.strength.value]
        new_comment = "A splicing assay was performed with allele_specific quantification using patient RNA showing <90% and >80% proportion of non-functional transcript."
        update_rule_result(result, strength=strength, comment=new_comment)
        return result
    if minigene_90 and not minigene:
        strength = strength_adjustment_by_one[result.strength.value]
        new_comment = "A mingene assay was performed showing >=90% proportion of non-functional transcript."
        update_rule_result(result, strength=strength, comment=new_comment)
        return result
    if minigene:
        strength = strength_adjustment_by_two[result.strength.value]
        new_comment = "A minigene assay was performed showing <90% and >80% proportion of non-functional transcript."
        update_rule_result(result, strength=strength, comment=new_comment)
        return result
    if patient_no_allele_specific_patho:
        strength = strength_adjustment_no_allelic_quant[result.strength.value]
        new_comment = "A splicing assay was performed with quantification (no allele-speicific quantification) using patient RNA showing >=50% proportion of non-functional transcript from both alleles."
        update_rule_result(result, strength=strength, comment=new_comment)
        return result
    return result


def update_rule_result(
    rule_result: RuleResult,
    comment: str,
    status: Optional[bool] = None,
    strength: Optional[evidence_strength] = None,
):
    """
    Update rule result based on classification
    """
    if status is not None:
        rule_result.status = status
    if strength is not None:
        rule_result.strength = strength
    rule_result.comment = rule_result.comment + " " + comment


def assess_splicing_data_bp7(
    rna_data: list[RNAData],
) -> tuple[bool, bool, str]:
    """
    Assess splicing data for applicability to BP7
    """
    performed = True
    no_quantification = 0
    bp7_allele, bp7_allele_not, bp7, bp7_no = False, False, False, False
    for entry in rna_data:
        if entry.quantification is None:
            no_quantification = +1
        elif (
            entry.patient_rna and entry.allelic.value is ALLELIC.TRUE.value
        ) or entry.minigene:
            if entry.quantification < 0.7:
                bp7_allele = True
            else:
                bp7_allele_not = True
        elif entry.patient_rna and entry.allelic.value is ALLELIC.FALSE.value:
            if entry.quantification <= 0.05:
                bp7 = True
            else:
                bp7_no = True
    if not rna_data:
        performed = False
        result = False
        comment = "No RNA assay performed."
    elif no_quantification == len(rna_data):
        result = False
        comment = "No quantification available for the RNA assay. Please check results of RNA assay manually."
    elif (bp7_allele_not or bp7_no) and (bp7 or bp7_allele):
        result = False
        comment = "Multiple RNA Assays were performed with contradictory or inconclusive results. Please check results of RNA assay manually."
    elif bp7_allele_not:
        result = False
        comment = "The performed RNA assays with allele specific quantification or minigene assay do not show the variant to have no effect on splicing."
    elif bp7_no:
        result = False
        comment = "The performed RNA assays without allele specific quantification do not show the variant to have no effect splicing."
    elif bp7_allele:
        result = True
        comment = "The performed RNA assay with allel specific quantification or minigene assay shows the variant not to have an effect on splicing."
    elif bp7:
        result = True
        comment = "The performed RNA assay without allel specific quantification shows the variant not to have an effect on splicing."
    else:
        result = False
        comment = (
            "RNA assay is inconclusive. Please check results of RNA assay manually."
        )
    return performed, result, comment
