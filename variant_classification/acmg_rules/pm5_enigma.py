#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info


class Pm5_ptc_enigma(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Implementing the rule specifications for ATM
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.PM5_RESULTS_PTC,),
        )

    @classmethod
    def assess_rule(cls, pm5_result: Optional[RuleResult]) -> RuleResult:
        if pm5_result:
            return pm5_result
        else:
            return RuleResult(
                "PM5",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                "PM5 not applicable to this variant.",
            )
