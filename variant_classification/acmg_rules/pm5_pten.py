#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from clinvar_utils import ClinVar, ClinVar_Status
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds


class Pm5_protein_pten(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Including the gene specifications for PTEN
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN_SIMILARITY,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_diff_aa: ClinVar,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if (
            clinvar_diff_aa.pathogenic
            and clinvar_diff_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            if num_thresholds_met is None:
                result = True
                comment = f"ATTENTION: No splicing prediction is available for variant under assessment. The following ClinVar entries show an amino acid change in the same position as pathogenic: {clinvar_diff_aa.ids}."
            elif num_thresholds_met > 0:
                result = True
                comment = f"The following ClinVar entries show an amino acid change in the same position as pathogenic: {clinvar_diff_aa.ids}."
            else:
                result = False
                comment = f"Variant is not predicted to not affect splicing. PM5_protein is therefore not applicable."
        else:
            result = False
            comment = "No ClinVar entries found that show an amino acid change in the same position as pathogenic."
        return RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
