#!/usr/bin/env python3
#
from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds


class Pm5_protein(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_results: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        clinvar_diff_aa = clinvar_results[ClinVar_Type.DIFF_AA_CHANGE]
        if clinvar_diff_aa.pathogenic:
            if num_thresholds_met is None:
                result = True
                comment = f"ATTENTION: No splicing prediction is available for variant under assessment. The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_diff_aa.ids}."
            elif num_thresholds_met > 0:
                result = True
                comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_diff_aa.ids}."
            else:
                result = False
                comment = f"Variant is not predicted to not affect splicing. PS1_protein is therefore not applicable."
            if clinvar_diff_aa.associated_ids and result:
                comment = (
                    comment
                    + f" The following ClinVar entries show the same amino acid exchange as likely pathogenic: {clinvar_diff_aa.associated_ids}."
                )
        else:
            result = False
            comment = "No ClinVar entries found that show an amino acid change in the same position as (likely) pathogenic."
        return RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm5_protein_pathogenic(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Assert that ClinVar variant is pathogenic, likely pathogenic would not be sufficient
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_results: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        clinvar_diff_aa = clinvar_results[ClinVar_Type.DIFF_AA_CHANGE]
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if (
            clinvar_diff_aa.pathogenic
            and clinvar_diff_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            if num_thresholds_met is None:
                result = True
                comment = f"ATTENTION: No splicing prediction is available for variant under assessment. The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_diff_aa.ids}."
            elif num_thresholds_met > 0:
                result = True
                comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_diff_aa.ids}."
            else:
                result = False
                comment = f"Variant is not predicted to not affect splicing. PS1_protein is therefore not applicable."
            if clinvar_diff_aa.associated_ids and result:
                comment = (
                    comment
                    + f" The following ClinVar entries show the same amino acid exchange as likely pathogenic: {clinvar_diff_aa.associated_ids}."
                )
        else:
            comment = "No ClinVar entries found that show an amino acid change in the same position as pathogenic."
            result = False
        return RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
