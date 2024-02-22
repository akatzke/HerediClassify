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


class Pp4(abstract_rule):
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
        if multifactorial_likelihood.multifactorial_likelihood is None:
            result = False
            comment = "No multifactorial likelihood given for variant."
            strength = evidence_strength.SUPPORTING

        num_thresholds_met = assess_thresholds(
            thresholds, multifactorial_likelihood.co_segregation
        )
        if num_thresholds_met is None:
            result = False
            comment = "No multifactorial likelihood given for variant."
            strength = evidence_strength.SUPPORTING

        elif num_thresholds_met == 0:
            result = False
            strength = evidence_strength.SUPPORTING
            comment = f"Multifactorial likelihood of {multifactorial_likelihood.multifactorial_likelihood} given for the variant meets no threshold for pathogenic evidence."
        else:
            result = True
            strength = thresholds.strengths[num_thresholds_met - 1]
            comment = f"Multifactorial likelihood of {multifactorial_likelihood.multifactorial_likelihood} given for variant meets threshold for {strength.value} pathogenic evidence."

        return RuleResult(
            "PP4",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
