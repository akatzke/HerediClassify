#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from acmg_rules.functional_splicing_assay_utils import (
    summarise_func_data,
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
        func_data: list[FunctionalData],
    ) -> RuleResult:
        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)
        if not func_data:
            result = False
            comment = f"No functional assay performed."
        elif pathogenic_count == 0:
            result = False
            comment = f"None of the {len(func_data)} functional assay(s) performed indicate pathogenicity of the variant."
        elif benign_count > 0 and pathogenic_count > 0:
            result = False
            comment = f"{pathogenic_count} of the {len(func_data)} perforemd assays indicate that the variant is pathogenic and {benign_count} of the {len(func_data)} indicate that the variant is benign. Due to this conflicting evidenced PS3 can not be applied."
        elif pathogenic_count > 0:
            result = True
            comment = f"{benign_count} of the {len(func_data)} performed assays indicate that the variant is pathogenic."
            if uncertain_count != 0:
                comment = (
                    comment
                    + f" ATTENTION: {uncertain_count} of the {len(func_data)} performed assays show no clear result."
                )
        else:
            result = False
            comment = (
                "Functional assay performed but results do not indicate pathogenicity."
            )
        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Ps3_prot_and_splice_assay(abstract_rule):
    """
    PS3: Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
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
            if func_data.pathogenic and splice_data.pathogenic:
                result = True
                comment = (
                    comment
                    + "Both the splice assay and the functional assay indicate pathogenicity."
                )
                type = rule_type.GENERAL
            elif func_data.pathogenic:
                result = True
                comment = comment + " The functional assay indicates pathogenicity."
                type = rule_type.PROTEIN
            elif splice_data.pathogenic:
                result = True
                comment = comment + " The splice assay indicates pathogenicity."
                type = rule_type.SPLICING
            else:
                result = False
                comment = (
                    comment
                    + " Neither functional assay nor splice assay indicate pathogenicity."
                )
                type = rule_type.GENERAL
        elif splice_data.performed:
            type = rule_type.SPLICING
            if splice_data.pathogenic:
                result = True
                comment = f"Splice assay performed and result indicates pathogenicity."
            else:
                result = False
                comment = f"Splice assay performed and result does not indicate pathogenicity."
        elif func_data.performed:
            type = rule_type.PROTEIN
            if func_data.pathogenic:
                result = True
                comment = (
                    f"Functional assay performed and result indicates pathogenicity."
                )
            else:
                result = False
                comment = f"Functional assay performed and result does not indicate pathogenicity."
        else:
            type = rule_type.GENERAL
            result = False
            comment = (
                "Functional assay performed but results do not indicate pathogenicity."
            )
        return RuleResult(
            "BS3",
            type,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Ps3_only_splice(abstract_rule):
    """
    PS3: Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
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
        elif splice_data.pathogenic:
            result = True
            comment = f"Splice assay performed and result indicates pathogenicity."
        else:
            result = False
            comment = (
                "Splice assay performed but results do not indicate pathogenicity."
            )
        return RuleResult(
            "BS3",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )
