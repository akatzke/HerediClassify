#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds


class Pp4_pcd(abstract_rule):
    """
    PP4: Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology.
    Here used for multifactorial likelihood
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_PICADAR_SCORE,
                class_info.VARIANT_MODIFIED_PICADAR_SCORE,
                class_info.THRESHOLD_PICADAR_SCORE,
                class_info.THRESHOLD_MODIFIED_PICADAR_SCORE,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        picadar_score: Optional[int],
        modified_picadar_score: Optional[int],
        threshold_picadar_score: Threshold,
        threshold_modified_picadar_score: Threshold,
    ) -> RuleResult:
        if isinstance(picadar_score,int) and isinstance(modified_picadar_score, int):
            num_thresholds_met = assess_thresholds(threshold_picadar_score, picadar_score)
            if num_thresholds_met == 0:
                result = False
                strength = evidence_strength.SUPPORTING
                comment = f"Patient has PICADAR score of {picadar_score}, therefore PP4 does not apply."
            else:
                result = True
                strength = threshold_picadar_score.strengths[num_thresholds_met - 1]
                comment = f"Patient has PICADAR score of {picadar_score}, therefore PP4 applies at evidence strength {strength.value}."
        elif isinstance(picadar_score,int):
            num_thresholds_met = assess_thresholds(threshold_picadar_score, picadar_score)
            if num_thresholds_met == 0:
                result = False
                strength = evidence_strength.SUPPORTING
                comment = f"Patient has PICADAR score of {picadar_score}, therefore PP4 does not apply."
            else:
                result = True
                strength = threshold_picadar_score.strengths[num_thresholds_met - 1]
                comment = f"Patient has PICADAR score of {picadar_score}, therefore PP4 applies at evidence strength {strength.value}."
        elif isinstance(modified_picadar_score,int):
            num_thresholds_met = assess_thresholds(threshold_modified_picadar_score, modified_picadar_score)
            if num_thresholds_met == 0:
                result = False
                strength = evidence_strength.SUPPORTING
                comment = f"Patient has PICADAR score of {modified_picadar_score}, therefore PP4 does not apply."
            else:
                result = True
                strength = threshold_picadar_score.strengths[num_thresholds_met - 1]
                comment = f"Patient has PICADAR score of {modified_picadar_score}, therefore PP4 applies at evidence strength {strength.value}."
        else:
            comment = f"Neither PICADAR score, nor modified PICADAR score were provided. PP4 does not apply."
            strength = evidence_strength.SUPPORTING
            result = False
        return RuleResult(
            "PP4",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
