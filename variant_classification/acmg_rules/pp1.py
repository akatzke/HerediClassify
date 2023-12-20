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
from acmg_rules.computation_evidence_utils import Threshold_evidence_strength
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
        thresholds: Threshold_evidence_strength,
    ) -> RuleResult:
        if multifactorial_likelihood.co_segregation is None:
            result = False
            comment = "No Co-segregation given for variant."
            strength = evidence_strength.SUPPORTING
        elif (
            thresholds.threshold_very_strong is not None
            and multifactorial_likelihood.co_segregation
            > thresholds.threshold_very_strong
        ):
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Co-segregation given for variant {multifactorial_likelihood.co_segregation} meets threshold for very strong pathogenic evidence strength {thresholds.threshold_very_strong}."
        elif (
            thresholds.threshold_strong is not None
            and multifactorial_likelihood.co_segregation > thresholds.threshold_strong
        ):
            result = True
            strength = evidence_strength.STRONG
            comment = f"Co-segregation given for variant {multifactorial_likelihood.co_segregation} meets threshold for strong pathogenic evidence strength {thresholds.threshold_strong}."
        elif (
            thresholds.threshold_moderate is not None
            and multifactorial_likelihood.co_segregation > thresholds.threshold_moderate
        ):
            result = True
            strength = evidence_strength.MODERATE
            comment = f"Co-segregation given for variant {multifactorial_likelihood.co_segregation} meets threshold for moderate pathogenic evidence strength {thresholds.threshold_moderate}."
        elif (
            thresholds.threshold_supporting is not None
            and multifactorial_likelihood.co_segregation
            > thresholds.threshold_supporting
        ):
            result = True
            strength = evidence_strength.SUPPORTING
            comment = f"Co-segregation given for variant {multifactorial_likelihood.co_segregation} meets threshold for supporting pathogenic evidence strength {thresholds.threshold_supporting}."
        else:
            result = False
            strength = evidence_strength.SUPPORTING
            comment = f"Co-segregation given for variant {multifactorial_likelihood.co_segregation} does not meet any threshold for pathogenic evidence."

        return RuleResult(
            "PP1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
