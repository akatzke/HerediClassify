#!/usr/bin/env python3
#
from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar


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
            comment = f"The following ClinVar entries show an amino acid change in the same position as (likely) pathogenic: {clinvar_diff_aa.pathogenic}."
            result = True
        else:
            comment = "No ClinVar entries found that show an amino acid change in the same position as (likely) pathogenic."
            result = False
        return RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm5_protein_pathogenic(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Assert that ClinVar variant is pathogenic, likely pathogenic would not be sufficient
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
        if (
            clinvar_diff_aa.pathogenic
            and clinvar_diff_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show an amino acid change in the same position as pathogenic: {clinvar_diff_aa.pathogenic}."
            result = True
        else:
            comment = "No ClinVar entries found that show an amino acid change in the same position as pathogenic."
            result = False
        return RuleResult(
            "PM5",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )
