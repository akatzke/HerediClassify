#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
)
from acmg_rules.pvs1 import Pvs1
from acmg_rules.pvs1_brca2 import Pvs1_brca2
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)


class Pvs1_atm(Pvs1):
    """
    PVS1: loss of function
    Following VCEP guidelines for ATM
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.ANNOTATED_TRANSCRIPT_LIST,
                class_info.VARIANT,
                class_info.POS_LAST_KNOWN_PATHO_PTC,
                class_info.SPLICE_RESULT,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript: list[TranscriptInfo],
        variant: VariantInfo,
        pos_last_known_patho_ptc: int,
        splice_result: Optional[RuleResult],
    ):
        if len(annotated_transcript) != 1:
            raise ValueError(
                "For ATM more than one transcript is being used for assessment of PVS1, despite only one disease relevant transcript being defined."
            )
        transcript = annotated_transcript[0]
        if type(transcript) is TranscriptInfo_exonic:
            result = Pvs1_brca2.assess_pvs1_frameshift_PTC_brca2(
                transcript, pos_last_known_patho_ptc
            )
        elif type(transcript) is TranscriptInfo_intronic:
            if splice_result is None:
                result = cls.assess_pvs1_splice(transcript)
            else:
                result = splice_result
        elif type(transcript) is TranscriptInfo_start_loss:
            result = cls.assess_pvs1_start_loss_atm(transcript)
        else:
            comment = f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            result = RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                comment,
            )
        return result

    @classmethod
    def assess_pvs1_start_loss_atm(
        cls, transcript: TranscriptInfo_start_loss
    ) -> RuleResult:
        """
        Assess PVS1 for start lost variants
        """
        if transcript.is_truncated_region_disease_relevant:
            comment = f"Alternative start codon leads to the exclusion of a disease relevant region."
            comment = comment + " " + transcript.comment_truncated_region
            result = True
            strength = evidence_strength.VERY_STRONG
        else:
            comment = f"Alternative start codon does not lead to the exclusion of a disease relevant region."
            result = False
            strength = evidence_strength.VERY_STRONG
        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
