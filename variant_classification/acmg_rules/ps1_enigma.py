#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Classification_Info, Info
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar


class Ps1_protein_enigma(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        clinvar_same_aa = clinvar_result[ClinVar_Type.SAME_AA_CHANGE]
        if num_thresholds_met is None:
            result = False
            strength = evidence_strength.STRONG
            comment = "No splicing prediction is available. Therefore PS1_protein can not be evaluated."
        elif (
            clinvar_same_aa.pathogenic
            and num_thresholds_met > 0
            and clinvar_same_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}."
            strength = evidence_strength.STRONG
            result = True
            if clinvar_same_aa.associated_ids:
                comment = (
                    comment
                    + f" The following ClinVar entries show the same amino acid change as likely pathogenic: {clinvar_same_aa.associated_ids}."
                )
        elif (
            clinvar_same_aa.pathogenic
            and num_thresholds_met > 0
            and clinvar_same_aa.highest_classification
            == ClinVar_Status.LIKELY_PATHOGENIC
        ):
            comment = f"The following ClinVar entries show the same amino acid change as likely pathogenic: {clinvar_same_aa.ids}."
            strength = evidence_strength.MODERATE
            result = True
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to not affect splicing."
            strength = evidence_strength.STRONG
            result = False
        else:
            comment = "No ClinVar entries found that show the same amino acid change as pathogenic."
            strength = evidence_strength.STRONG
            result = False
        return RuleResult(
            "PS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
