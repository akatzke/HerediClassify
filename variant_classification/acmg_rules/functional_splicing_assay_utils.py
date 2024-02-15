#!/usr/bin/env python3

import logging

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
    minigene_90, minigene, patient_90, patient, patient_no_allele_specific = (
        0,
        0,
        0,
        0,
        0,
    )
    for entry in rna_data:
        if not entry.quantification or entry.quantification is None:
            no_quantification += 1
        if entry.minigene and entry.allelic is ALLELIC.CONSTRUCT:
            if entry.quantification >= 0.9:
                minigene_90 += 1
            elif entry.quantification < 0.9 and entry.quantification > 0.8:
                minigene += 0
        elif entry.patient_rna:
            if entry.allelic is ALLELIC.TRUE:
                if entry.quantification >= 0.9:
                    patient_90 += 1
                elif entry.quantification < 0.9 and entry.quantification > 0.8:
                    patient += 1
            elif entry.allelic is ALLELIC.FALSE and entry.quantification >= 0.5:
                patient_no_allele_specific += 1
    strength_adjustment_by_one = {
        evidence_strength.VERY_STRONG: evidence_strength.STRONG,
        evidence_strength.STRONG: evidence_strength.MODERATE,
        evidence_strength.MODERATE: evidence_strength.SUPPORTING,
        evidence_strength.SUPPORTING: evidence_strength.SUPPORTING,
    }
    strength_adjustment_by_two = {
        evidence_strength.VERY_STRONG: evidence_strength.MODERATE,
        evidence_strength.STRONG: evidence_strength.SUPPORTING,
        evidence_strength.MODERATE: evidence_strength.SUPPORTING,
        evidence_strength.SUPPORTING: evidence_strength.SUPPORTING,
    }
    strength_adjustment_no_allelic_quant = {
        evidence_strength.VERY_STRONG: evidence_strength.STRONG,
        evidence_strength.STRONG: evidence_strength.MODERATE,
        evidence_strength.MODERATE: evidence_strength.MODERATE,
        evidence_strength.SUPPORTING: evidence_strength.MODERATE,
    }
    if not rna_data:
        result.comment = result.comment + " No RNA assay performed."
    elif no_quantification == len(rna_data):
        result.comment = (
            result.comment
            + " No quantification available for the RNA assay. Please check resulst of RNA assay manually."
        )
    if patient_90 > 0 and patient == 0:
        result.name = "PVS1_RNA"
        result.comment = (
            result.comment
            + " A splicing assay was performed with allele_specific quantification using patient RNA showing >=90% proportion of non-functional transcript."
        )
    elif patient > 0:
        result.name = "PVS1_RNA"
        result.strength = strength_adjustment_by_one[result.strength]
        result.comment = (
            result.comment
            + " A splicing assay was performed with allele_specific quantification using patient RNA showing <90% and >80% proportion of non-functional transcript."
        )
    elif minigene_90 > 0 and minigene == 0:
        result.name = "PVS1_RNA"
        result.strength = strength_adjustment_by_one[result.strength]
        result.comment = (
            result.comment
            + " A mingene assay was performed showing >=90% proportion of non-functional transcript."
        )
    elif minigene:
        result.name = "PVS1_RNA"
        result.strength = strength_adjustment_by_two[result.strength]
        result.comment = (
            result.comment
            + " A minigene assay was performed showing <90% and >80% proportion of non-functional transcript."
        )
    elif patient_no_allele_specific > 0:
        result.strength = strength_adjustment_no_allelic_quant[result.strength]
        result.comment = (
            result.comment
            + " A splicing assay was performed with quantification (no allele-speicific quantification) using patient RNA showing >=50% proportion of non-functional transcript from both alleles."
        )
    return result


def assess_splicing_data_bp7(
    rna_data: list[RNAData],
) -> tuple[bool, bool, str]:
    """
    Assess splicing data for applicability to BP7
    """
    performed = True
    no_quantification = 0, 0
    bp7_allele, bp7_allele_not, bp7, bp7_no = 0, 0, 0, 0
    for entry in rna_data:
        if not entry.quantification or entry.quantification is None:
            no_quantification = +1
        elif (entry.patient_rna and entry.allelic is ALLELIC.TRUE) or entry.minigene:
            if entry.quantification < 0.7:
                bp7_allele += 1
            else:
                bp7_allele_not += 1
        elif entry.patient_rna and entry.allelic is ALLELIC.FALSE:
            if entry.quantification == 0:
                bp7 += 1
            else:
                bp7_no += 1
    if not rna_data:
        performed = False
        result = False
        comment = "No RNA assay performed."
    elif no_quantification == len(rna_data):
        result = False
        comment = "No quantification available for the RNA assay. Please check results of RNA assay manually."
    elif (bp7_allele_not > 0 or bp7_no > 0) and (bp7 > 0 or bp7_allele > 0):
        result = False
        comment = "Multiple RNA Assays were performed with contradictory results. Please check results manually."
    elif bp7_allele_not > 0:
        result = False
        comment = "The performed RNA assays with allele specific quantification or minigene assay, do not show the variant to have no effect on splicing."
    elif bp7_no > 0:
        result = False
        comment = "The performed RNA assays without allele specific quantification, do not show the variant to have no effect splicing."
    elif bp7_allele:
        result = True
        comment = "The performed RNA assay with allel specific quantification or minigene assay, shows the variant not to have an effect on splicing."
    elif bp7:
        result = True
        comment = "The performed RNA assay without allel specific quantification, shows the variant not to have an effect on splicing."
    else:
        result = False
        comment = "RNA assay is inconclusive. Please check results manually."
    return performed, result, comment
