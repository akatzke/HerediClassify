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
from acmg_rules.functional_splicing_assay_utils import (
    summarise_func_data,
)


class Bs3(abstract_rule):
    """
    BS3: Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing.
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
        func_data: list[FunctionalData],
    ) -> RuleResult:
        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)
        if not func_data:
            result = False
            comment = f"No functional assay performed."
        elif benign_count == 0:
            result = False
            comment = f"None of the {len(func_data)} functional assay(s) performed indicate benignity of the variant."
        elif benign_count > 0 and pathogenic_count > 0:
            result = False
            comment = f"{pathogenic_count} of the {len(func_data)} perforemd assays indicate that the variant is pathogenic and {benign_count} of the {len(func_data)} indicate that the variant is benign. Due to this conflicting evidenced BS3 can not be applied."
        elif benign_count > 0:
            result = True
            comment = f"{benign_count} of the {len(func_data)} performed assays indicate that the variant is benign."
            if uncertain_count != 0:
                comment = (
                    comment
                    + f" ATTENTION: {uncertain_count} of the {len(func_data)} performed assays show no clear result."
                )
        else:
            result = False
            comment = (
                "Functional assay performed but results do not indicate benignity."
            )
        return RuleResult(
            "BS3",
            rule_type.PROTEIN,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )
