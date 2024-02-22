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


class Bs4(abstract_rule):
    """
    BS4: Lack of segregation in affected members of a family.
    Here used for Co-segregation and multifactorial likelihood
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
        if (
            multifactorial_likelihood.co_segregation is None
            and multifactorial_likelihood.multifactorial_likelihood is None
        ):
            return RuleResult(
                "BS4",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                False,
                evidence_strength.STRONG,
                "No Co-segregation or multifactorial likelihood given for variant.",
            )
        elif multifactorial_likelihood.multifactorial_likelihood is not None:
            likelihood = multifactorial_likelihood.multifactorial_likelihood
            name = "Multifactorial likelihood"
        elif multifactorial_likelihood.co_segregation is not None:
            likelihood = multifactorial_likelihood.co_segregation
            name = "Co-segregation"
        else:
            return RuleResult(
                "BS4",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                False,
                evidence_strength.STRONG,
                "No Co-segregation or multifactorial likelihood given for variant.",
            )
        num_thresholds_met = assess_thresholds(threshold, likelihood)
        if num_thresholds_met is None:
            result = False
            strength = evidence_strength.STRONG
            comment = f"No {name} was provided. BS4 could not be evaluated."
        elif num_thresholds_met == 0:
            result = False
            strength = evidence_strength.STRONG
            comment = f"{name} of {likelihood} given for variant does not meet any threshold for benign evidence."
        else:
            result = True
            strength = threshold.strengths[num_thresholds_met - 1]
            comment = f"{name} of {likelihood} given for variant meets threshold for {strength.value} benign evidence."
        return RuleResult(
            "BS4",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )
