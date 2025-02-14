#!/usr/bin/env python3

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


class Pp1(abstract_rule):
    """
    PP1: Co-segregation with disease in multiple affected family members in a gene definitively known to cause the disease.
    Here used for Co-segregation
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_MULTIFACTORIAL_LIKELIHOOD,
                class_info.THRESHOLD_LIKELIHOOD_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        multifactorial_likelihood: MultifactorialLikelihood,
        thresholds: Threshold,
    ) -> RuleResult:
        num_thresholds_met = assess_thresholds(
            thresholds, multifactorial_likelihood.co_segregation
        )
        if num_thresholds_met is None:
            result = False
            comment = "No co-segregation data available for variant."
            strength = evidence_strength.SUPPORTING

        elif num_thresholds_met == 0:
            result = False
            strength = evidence_strength.SUPPORTING
            comment = f"Co-segregation of {multifactorial_likelihood.co_segregation} given for the variant meets no threshold for pathogenic evidence."
        else:
            result = True
            strength = thresholds.strengths[num_thresholds_met - 1]
            comment = f"Co-segregation of {multifactorial_likelihood.co_segregation} given for variant meets threshold for {strength.value} pathogenic evidence ({thresholds.thresholds[num_thresholds_met -1]})."
        return RuleResult(
            "PP1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
