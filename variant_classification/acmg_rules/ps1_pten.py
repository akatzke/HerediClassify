#!/usr/bin/env python3


from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Classification_Info, Info
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar


class Ps1_splicing_pten(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR_SPLICEAI_SPLICE,),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
    ) -> RuleResult:
        clinvar_same_nucleotide = clinvar_result[ClinVar_Type.SAME_NUCLEOTIDE]
        clinvar_same_splice_site = clinvar_result[ClinVar_Type.SAME_SPLICE_SITE]
        if (
            clinvar_same_nucleotide.pathogenic
            and clinvar_same_nucleotide.highest_classification
            == ClinVar_Status.PATHOGENIC
        ) or (
            clinvar_same_splice_site.pathogenic
            and clinvar_same_splice_site.highest_classification
            == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show splice variants at the same splice site to be pathogenic: {clinvar_same_nucleotide.ids}, {clinvar_same_splice_site.ids}."
            result = True
        else:
            comment = "No ClinVar entries found that show splice variants at the same nucleotide position as pathogenic."
            result = False
        return RuleResult(
            "PS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )
