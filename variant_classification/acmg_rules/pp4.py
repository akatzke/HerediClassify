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
            (class_info.VARIANT_MULTIFACTORIAL_LIKELIHOOD,),
        )

    @classmethod
    def assess_rule(
        cls,
        multifactorial_likelihood: MultifactorialLikelihood,
        thresholds: Threshold_evidence_strength,
    ) -> RuleResult:
        if multifactorial_likelihood.multifactorial_likelihood is None:
            result = False
            comment = "No multifactorial likelihood given for variant."
            strength = evidence_strength.SUPPORTING
        elif (
            thresholds.threshold_very_strong is not None
            and multifactorial_likelihood.multifactorial_likelihood
            > thresholds.threshold_very_strong
        ):
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Multifactorial likelihood given for variant {multifactorial_likelihood.multifactorial_likelihood} meets threshold for very strong pathogenic evidence strength {thresholds.threshold_very_strong}."
        elif (
            thresholds.threshold_strong is not None
            and multifactorial_likelihood.multifactorial_likelihood
            > thresholds.threshold_strong
        ):
            result = True
            strength = evidence_strength.STRONG
            comment = f"Multifactorial likelihood given for variant {multifactorial_likelihood.multifactorial_likelihood} meets threshold for strong pathogenic evidence strength {thresholds.threshold_strong}."
        elif (
            thresholds.threshold_moderate is not None
            and multifactorial_likelihood.multifactorial_likelihood
            > thresholds.threshold_moderate
        ):
            result = True
            strength = evidence_strength.MODERATE
            comment = f"Multifactorial likelihood given for variant {multifactorial_likelihood.multifactorial_likelihood} meets threshold for moderate pathogenic evidence strength {thresholds.threshold_moderate}."
        elif (
            thresholds.threshold_supporting is not None
            and multifactorial_likelihood.multifactorial_likelihood
            > thresholds.threshold_supporting
        ):
            result = True
            strength = evidence_strength.SUPPORTING
            comment = f"Multifactorial likelihood given for variant {multifactorial_likelihood.multifactorial_likelihood} meets threshold for supporting pathogenic evidence strength {thresholds.threshold_supporting}."
        else:
            result = False
            strength = evidence_strength.SUPPORTING
            comment = f"Multifactorial likelihood given for variant {multifactorial_likelihood.multifactorial_likelihood} does not meet any threshold for pathogenic evidence."

        return RuleResult(
            "PP4",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
