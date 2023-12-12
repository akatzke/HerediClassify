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
from variant_classification.var_type import VARTYPE


class Pvs1_cdh1(Pvs1):
    """
    PVS1: loss of function
    Following VCEP guidelines for CDH1
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
                class_info.SPLICE_RESULT,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript: list[TranscriptInfo],
        variant: VariantInfo,
        splice_result: Optional[RuleResult],
    ):
        if len(annotated_transcript) != 1:
            raise ValueError(
                "For CDH1 more than one transcript is being used for assessment of PVS1, despite only one disease relevant transcript being defined."
            )
        transcript = annotated_transcript[0]
        if type(transcript) is TranscriptInfo_exonic:
            result = cls.assess_pvs1_frameshift_PTC_cdh1(transcript)
        elif type(transcript) is TranscriptInfo_intronic:
            if splice_result is None:
                result = cls.assess_pvs1_splice(transcript)
            else:
                result = splice_result
        elif type(transcript) is TranscriptInfo_start_loss:
            result = cls.assess_pvs1_start_loss(transcript)
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
    def assess_pvs1_frameshift_PTC_cdh1(
        cls, transcript: TranscriptInfo_exonic
    ) -> RuleResult:
        if VARTYPE.STOP_GAINED in transcript.var_type:
            if transcript.var_start >= 4 and transcript.var_start <= 2388:
                comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
                result = True
                strength = evidence_strength.VERY_STRONG
            elif transcript.var_start >= 2389 and transcript.var_start <= 2508:
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Truncated region is disease relevant."
                result = True
                strength = evidence_strength.STRONG
            elif transcript.var_start >= 2509 and transcript.var_start <= 2646:
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Truncated region is not disease relevant and less than 10% of protein length are lost."
                result = True
                strength = evidence_strength.MODERATE
            else:
                raise ValueError(
                    f"The start position of the variant {transcript.var_start} in transcript {transcript.transcript_id} lies outside the range for nonsense and framshift variants of c.4-c.2646."
                )
        elif transcript.frameshift == -1:
            if transcript.var_start >= 4 and transcript.var_start <= 2375:
                comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
                result = True
                strength = evidence_strength.VERY_STRONG
            elif transcript.var_start >= 2376 and transcript.var_start <= 2443:
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Truncated region is disease relevant."
                result = True
                strength = evidence_strength.STRONG
            elif transcript.var_start >= 2444 and transcript.var_start <= 2646:
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Truncated region is not disease relevant and less than 10% of protein length are lost."
                result = True
                strength = evidence_strength.MODERATE
            else:
                raise ValueError(
                    f"The start position of the variant {transcript.var_start} in transcript {transcript.transcript_id} lies outside the range for nonsense and framshift variants of c.4-c.2646."
                )
        elif transcript.frameshift == +1:
            if transcript.var_start >= 4 and transcript.var_start <= 2352:
                comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
                result = True
                strength = evidence_strength.VERY_STRONG
            elif transcript.var_start >= 2353 and transcript.var_start <= 2499:
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Truncated region is disease relevant."
                result = True
                strength = evidence_strength.STRONG
            elif transcript.var_start >= 2500 and transcript.var_start <= 2646:
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Truncated region is not disease relevant and less than 10% of protein length are lost."
                result = True
                strength = evidence_strength.MODERATE
            else:
                raise ValueError(
                    f"The start position of the variant {transcript.var_start} in transcript {transcript.transcript_id} lies outside the range for nonsense and framshift variants of c.4-c.2646."
                )
        else:
            comment = f"Variant in transcript {transcript.transcript_id} does not meet any of the specified criteria for PALB2."
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
