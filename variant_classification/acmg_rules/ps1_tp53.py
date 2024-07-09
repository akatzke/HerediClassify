#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from acmg_rules.functional_splicing_assay_utils import (
    assess_splicing_data_bp7,
)
from information import Classification_Info, Info
from clinvar_utils import ClinVar_Status, ClinVar_Type, ClinVar
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from variant import RNAData


class Ps1_protein_tp53(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splicing_assay: Optional[list[RNAData]],
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        clinvar_same_aa = clinvar_result[ClinVar_Type.SAME_AA_CHANGE]
        if splicing_assay:
            performed, result_assay, comment_assay = assess_splicing_data_bp7(
                splicing_assay
            )
        else:
            performed = False
            result_assay = False
        if (
            clinvar_same_aa.pathogenic
            and performed
            and result_assay
            and clinvar_same_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}. A splice assays shows that the variant does not affect splicing."
            strength = evidence_strength.STRONG
            result = True
            if clinvar_same_aa.associated_ids:
                comment = (
                    comment
                    + f" The following ClinVar entries show the same amino acid change as likely pathogenic: {clinvar_same_aa.associated_ids}."
                )
        elif (
            clinvar_same_aa.pathogenic
            and performed
            and not result_assay
            and clinvar_same_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"A splice assay shows that the variant affects splicing."
            strength = evidence_strength.STRONG
            result = False
        elif num_thresholds_met is None:
            result = False
            strength = evidence_strength.STRONG
            comment = "No splicing prediction is available. Therefore, PS1_protein can not be evaluated."
        elif (
            clinvar_same_aa.pathogenic
            and num_thresholds_met == 0
            and clinvar_same_aa.highest_classification == ClinVar_Status.PATHOGENIC
        ):
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}. SpliceAI does not predict an effect on splicing for this variant."
            strength = evidence_strength.MODERATE
            result = True
        elif num_thresholds_met > 0:
            comment = f"Variant is predicted to affect splicing."
            strength = evidence_strength.STRONG
            result = False
        else:
            comment = "No ClinVar entries found that show the same amino acid change as pathogenic."
            strength = evidence_strength.STRONG
            result = False
        return RuleResult(
            "PS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


class Ps1_splicing_tp53(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_SPLICEAI_SPLICE,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_result: dict[ClinVar_Type, ClinVar],
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        clinvar_same_nucleotide = clinvar_result[ClinVar_Type.SAME_NUCLEOTIDE]
        if num_thresholds_met is None:
            comment = "No splicing prediction is available. Therefore PS1_splicing can not be evaluated."
            result = False
        elif clinvar_same_nucleotide.pathogenic and num_thresholds_met > 0:
            comment = f"The following ClinVar entries show splice variants at the same nucleotide position to be pathogenic: {clinvar_same_nucleotide.ids}."
            result = True
        elif num_thresholds_met == 0:
            comment = f"The variant is not predicted to affect splicing."
            result = False
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
