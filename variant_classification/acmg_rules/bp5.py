#!/usr/bin/env python

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from variant import MultifactorialLikelihood


class Bp5_enigma(abstract_rule):
    """
    BP5: Variant found in a case with an alternative molecular basis for disease
    Modified for multifactorial likelihood data by ENIGMA Consortium
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_MULTIFACTORIAL_LIKELIHOOD,
                class_info.THRESHOLD_LIKELIHOOD_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        multifactorial_likelihood: MultifactorialLikelihood,
        threshold: Threshold,
    ) -> RuleResult:
        if multifactorial_likelihood.multifactorial_likelihood is None:
            result = False
            comment = "No multifactorial likelihood given for variant."
            strength = evidence_strength.SUPPORTING

        num_thresholds_met = assess_thresholds(
            threshold, multifactorial_likelihood.co_segregation
        )
        if num_thresholds_met is None:
            result = False
            comment = "No multifactorial likelihood given for variant."
            strength = evidence_strength.SUPPORTING

        elif num_thresholds_met == 0:
            result = False
            strength = evidence_strength.SUPPORTING
            comment = f"Multifactorial likelihood of {multifactorial_likelihood.multifactorial_likelihood} given for the variant meets no threshold for benign evidence."
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"Multifactorial likelihood of {multifactorial_likelihood.multifactorial_likelihood} given for variant meets threshold for {strength.value} beningn evidence ({threshold.thresholds[num_thresholds_met -1]})."

        return RuleResult(
            "BP5",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )
