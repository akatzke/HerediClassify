#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Info, Classification_Info
from variant import FunctionalData


class Ps3(abstract_rule):
    """
    PS3: Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY,),
        )

    @classmethod
    def assess_rule(
        cls,
        func_data: FunctionalData,
    ) -> RuleResult:
        if not func_data.performed:
            result = False
            comment = f"No functional assay performed."
        elif func_data.performed and func_data.pathogenic:
            result = True
            comment = f"Functional assay performed and result indicates benignity."
        else:
            result = False
            comment = (
                "Functional assay performed but results do not indicate benignity."
            )
        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )
