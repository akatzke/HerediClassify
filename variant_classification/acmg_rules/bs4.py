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
        thresholds: Threshold_evidence_strength,
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
        if (
            thresholds.threshold_very_strong is not None
            and likelihood > thresholds.threshold_very_strong
        ):
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"{name} given for variant {likelihood} meets threshold for very strong pathogenic evidence strength {thresholds.threshold_very_strong}."
        elif (
            thresholds.threshold_strong is not None
            and likelihood > thresholds.threshold_strong
        ):
            result = True
            strength = evidence_strength.STRONG
            comment = f"{name} given for variant {likelihood} meets threshold for strong pathogenic evidence strength {thresholds.threshold_strong}."
        elif (
            thresholds.threshold_moderate is not None
            and likelihood > thresholds.threshold_moderate
        ):
            result = True
            strength = evidence_strength.MODERATE
            comment = f"{name} given for variant {likelihood} meets threshold for moderate pathogenic evidence strength {thresholds.threshold_moderate}."
        elif (
            thresholds.threshold_supporting is not None
            and likelihood > thresholds.threshold_supporting
        ):
            result = True
            strength = evidence_strength.SUPPORTING
            comment = f"{name} given for variant {likelihood} meets threshold for supporting pathogenic evidence strength {thresholds.threshold_supporting}."
        else:
            result = False
            strength = evidence_strength.STRONG
            comment = f"{name} given for variant {likelihood} does not meet any threshold for pathogenic evidence."

        return RuleResult(
            "BS4",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )
