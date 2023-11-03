#!/usr/bin/env python3
#
from typing import Callable

from refactoring.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
import refactoring.information as info
from refactoring.clinvar_utils import CLINVAR_TYPE, ClinVar


class pm5_protein(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (info.classification_information.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_results: dict[CLINVAR_TYPE, ClinVar]) -> RuleResult:
        clinvar_diff_aa = clinvar_results[CLINVAR_TYPE.DIFF_AA_CHANGE]
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


class pm5_splicing(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    """

    @classmethod
    def get_assess_rule(
        cls,
    ) -> tuple[Callable, tuple[info.classification_information, ...]]:
        return (
            cls.assess_rule,
            (info.classification_information.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_results: dict[CLINVAR_TYPE, ClinVar]) -> RuleResult:
        clinvar_same_splice_site = clinvar_results[CLINVAR_TYPE.SAME_SPLICE_SITE]
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
