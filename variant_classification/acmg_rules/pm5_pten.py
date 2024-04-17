#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Info, Classification_Info
from clinvar_utils import ClinVar, ClinVar_Status


class Pm5_protein_pten(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    Including the gene specifications for PTEN
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN_SIMILARITY,
                class_info.VARIANT_HOTSPOT_ANNOTATION,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_diff_aa: ClinVar,
        variant_in_hotspot: bool,
    ) -> RuleResult:
        if variant_in_hotspot:
            return RuleResult(
                "PM5",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                comment=f"Variant located in mutational hotspot. PM5 is not applicable.",
            )
        if (
            clinvar_diff_aa.pathogenic
            and clinvar_diff_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show an amino acid change in the same position as pathogenic: {clinvar_diff_aa.ids}."
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
