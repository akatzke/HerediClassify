#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
)
from acmg_rules.pvs1 import Pvs1
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)


class Pvs1_brca1(Pvs1):
    """
    PVS1: loss of function
    Following VCEP guidelines for BRCA1
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
                "For BRCA1 more than one transcript is being used for assessment of PVS1, despite only one disease relevant transcript being defined."
            )
        transcript = annotated_transcript[0]
        if type(transcript) is TranscriptInfo_exonic:
            result = cls.assess_pvs1_frameshift_PTC_brca1(
                transcript, pos_last_known_patho_ptc
            )
        elif type(transcript) is TranscriptInfo_intronic:
            if splice_result is None:
                result = cls.assess_pvs1_splice(transcript)
            else:
                result = splice_result
        elif type(transcript) is TranscriptInfo_start_loss:
            result = cls.assess_pvs1_start_loss_pathogenic_very_strong()
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
    def assess_pvs1_frameshift_PTC_brca1(
        cls, transcript: TranscriptInfo_exonic, pos_last_known_patho_ptc: int
    ) -> RuleResult:
        if transcript.is_NMD:
            comment = (
                f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            )
            if transcript.is_affected_exon_disease_relevant:
                comment = comment + "Truncated region is disease relevant."
                result = True
            else:
                comment = comment + "Truncated region is not disease relevant."
                result = False
            strength = evidence_strength.VERY_STRONG
        else:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD."
            if transcript.var_start <= pos_last_known_patho_ptc:
                comment = (
                    comment
                    + "Truncated/altered region is critical to protein function."
                )
                result = True
                strength = evidence_strength.STRONG
            else:
                comment = (
                    "Role of truncated/alterend region in protein function is unknown."
                )
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
