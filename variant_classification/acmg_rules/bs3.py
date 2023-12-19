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
        func_data: FunctionalData,
    ) -> RuleResult:
        if not func_data.performed:
            result = False
            comment = f"No functional assay performed."
        elif func_data.benign:
            result = True
            comment = f"Functional assay performed and result indicates benignity."
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


class Bs3_prot_and_splice_assay(abstract_rule):
    """
    BS3: Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY, class_info.SPLICING_ASSAY),
        )

    @classmethod
    def assess_rule(
        cls, func_data: FunctionalData, splice_data: FunctionalData
    ) -> RuleResult:
        if not func_data.performed and not splice_data.performed:
            result = False
            comment = f"No functional assay or splice assay performed."
            type = rule_type.GENERAL
        elif func_data.performed and splice_data.performed:
            comment = f"Both functional and splice assay were performed."
            if func_data.benign and splice_data.benign:
                result = True
                comment = (
                    comment
                    + "Both the splice assay and the functional assay indicate benignity."
                )
                type = rule_type.GENERAL
            elif func_data.benign:
                result = True
                comment = comment + " The functional assay indicates benignity."
                type = rule_type.PROTEIN
            elif splice_data.benign:
                result = True
                comment = comment + " The splice assay indicates benignity."
                type = rule_type.SPLICING
            else:
                result = False
                comment = (
                    comment
                    + " Neither functional assay nor splice assay indicate benignity."
                )
                type = rule_type.GENERAL
        elif splice_data.performed:
            type = rule_type.SPLICING
            if splice_data.benign:
                result = True
                comment = f"Splice assay performed and result indicates benignity."
            else:
                result = False
                comment = (
                    f"Splice assay performed and result does not indicate benignity."
                )
        elif func_data.performed:
            type = rule_type.PROTEIN
            if func_data.benign:
                result = True
                comment = f"Functional assay performed and result indicates benignity."
            else:
                result = False
                comment = f"Functional assay performed and result does not indicate benignity."
        else:
            type = rule_type.GENERAL
            result = False
            comment = (
                "Functional assay performed but results do not indicate benignity."
            )
        return RuleResult(
            "BS3",
            type,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Bs3_only_splice(abstract_rule):
    """
    BS3: Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.SPLICING_ASSAY,),
        )

    @classmethod
    def assess_rule(
        cls,
        splice_data: FunctionalData,
    ) -> RuleResult:
        if not splice_data.performed:
            result = False
            comment = f"No splice assay performed."
        elif splice_data.benign:
            result = True
            comment = f"Splice assay performed and result indicates benignity."
        else:
            result = False
            comment = "Splice assay performed but results do not indicate benignity."
        return RuleResult(
            "BS3",
            rule_type.SPLICING,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )
