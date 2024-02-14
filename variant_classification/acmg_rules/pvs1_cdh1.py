#!/usr/bin/env python3


from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from acmg_rules.pvs1 import Pvs1
from information import Classification_Info, Info
from variant import RNAData, TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from acmg_rules.computation_evidence_utils import (
    Threshold,
)
from acmg_rules.functional_splicing_assay_utils import (
    adjust_strength_according_to_rna_data_pvs1,
)
from var_type import VARTYPE


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
                class_info.SPLICING_ASSAY,
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript: list[TranscriptInfo],
        variant: VariantInfo,
        splice_result: Optional[RuleResult],
        splice_assay: Optional[list[RNAData]],
        threshold_diff_len_prot_percent: float,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ):
        results = []
        for transcript in annotated_transcript:
            if isinstance(transcript, TranscriptInfo_exonic):
                result = cls.assess_pvs1_frameshift_PTC_cdh1(transcript)
                results.append(result)
            elif isinstance(transcript, TranscriptInfo_intronic):
                if splice_result is None:
                    splice_result = cls.assess_pvs1_splice(
                        transcript,
                        prediction_dict,
                        threshold,
                        threshold_diff_len_prot_percent,
                    )
                if splice_assay:
                    splice_result = adjust_strength_according_to_rna_data_pvs1(
                        splice_assay, splice_result
                    )
                results.append(splice_result)
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result = cls.assess_pvs1_start_loss(transcript)
                results.append(result)
        if len(results) == 0:
            comment = f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            result = RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                comment,
            )
        else:
            result = summarise_results_per_transcript(results)
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
