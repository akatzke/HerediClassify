#!/usr/bin/env python3
#
from typing import Callable

from variant_classification.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from variant_classification.information import Info, Classification_Info
from variant_classification.clinvar_utils import ClinVar_Type, ClinVar


class Pm5_protein(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_results: dict[ClinVar_Type, ClinVar]) -> RuleResult:
        clinvar_diff_aa = clinvar_results[ClinVar_Type.DIFF_AA_CHANGE]
        if clinvar_diff_aa.pathogenic:
            comment = f"The following ClinVar entries show an amino acid change in the same position as pathogenic: {clinvar_diff_aa.pathogenic}."
            result = True
        else:
            comment = "No matches found for variant."
            result = False
        return RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm5_splicing(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_results: dict[ClinVar_Type, ClinVar]) -> RuleResult:
        clinvar_same_splice_site = clinvar_results[ClinVar_Type.SAME_SPLICE_SITE]
        if clinvar_same_splice_site.pathogenic:
            comment = f"The following ClinVar entries show the variants in the same splice site as pathogenic: {clinvar_same_splice_site.pathogenic}."
            result = True
        else:
            comment = "No matches found for variant."
            result = False
        return RuleResult(
            "PM5",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
