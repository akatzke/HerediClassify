#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_prediction_tool
from acmg_rules.pvs1 import Pvs1
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo, FunctionalData
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from var_type import VARTYPE
from acmg_rules.computation_evidence_utils import Threshold


class Pvs1_palb2(Pvs1):
    """
    PVS1: loss of function
    Following VCEP guidelines for PALB2
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
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
                class_info.SPLICE_RESULT,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcript: list[TranscriptInfo],
        variant: VariantInfo,
        pos_last_known_patho_ptc_dict: dict[str, int],
        threshold_diff_len_prot_percent: float,
        splice_result: Optional[RuleResult],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splice_assay: FunctionalData,
    ):
        results = []
        for transcript in annotated_transcript:
            if isinstance(transcript, TranscriptInfo_exonic):
                try:
                    pos_last_known_patho_ptc = pos_last_known_patho_ptc_dict[
                        transcript.transcript_id
                    ]
                except KeyError:
                    raise KeyError(
                        f"Transcript {transcript.transcript_id} not in disease relevant transcripts: {pos_last_known_patho_ptc_dict.keys()}. Transcript should have been filtered out earlier."
                    )
                result = cls.assess_pvs1_frameshift_PTC_palb2(
                    transcript, pos_last_known_patho_ptc
                )
                results.append(result)
            elif isinstance(transcript, TranscriptInfo_intronic):
                if splice_assay.performed:
                    if splice_assay.pathogenic:
                        result = True
                        comment = f"A splice assay was performed showing a detrimental effect on splicing by the variant."
                    else:
                        result = False
                        comment = f"A splice assay was performed showing no detrimental effect on splicing by the variant."
                    return RuleResult(
                        "PVS1",
                        rule_type.SPLICING,
                        evidence_type.PATHOGENIC,
                        result,
                        evidence_strength.VERY_STRONG,
                        comment,
                    )
                if splice_result is None:
                    result = cls.assess_pvs1_splice_palb2(
                        transcript,
                        prediction_dict,
                        threshold,
                        threshold_diff_len_prot_percent,
                    )
                else:
                    result = splice_result
                results.append(result)
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result = cls.assess_pvs1_start_loss_pathogenic_very_strong()
                results.append(result)
        if len(results) == 0:
            comment = f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}."
            final_result = RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                comment,
            )
        else:
            final_result = summarise_results_per_transcript(results)
        return final_result

    @classmethod
    def assess_pvs1_frameshift_PTC_palb2(
        cls, transcript: TranscriptInfo_exonic, pos_last_known_patho_ptc: int
    ) -> RuleResult:
        if transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD and in a disease relevant transcript."
            result = True
            strength = evidence_strength.VERY_STRONG
        elif VARTYPE.STOP_GAINED in transcript.var_type:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Variant type is nonsense."
            if transcript.ptc <= pos_last_known_patho_ptc:
                comment = (
                    comment + f" Truncated region located in disease relevant region."
                )
                result = True
                strength = evidence_strength.VERY_STRONG
            else:
                comment = (
                    comment
                    + f" Truncated region not located in disease relevnat region."
                )
                result = True
                strength = evidence_strength.MODERATE
        elif VARTYPE.FRAMESHIFT_VARIANT in transcript.var_type:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Variant type is frameshift."
            if (
                transcript.is_truncated_region_disease_relevant
                and transcript.ptc <= pos_last_known_patho_ptc
            ):
                comment = (
                    comment
                    + f" Truncated region located in disease relevant region. PTC is located upstream of p.His1184."
                )
                result = True
                strength = evidence_strength.STRONG
            elif transcript.var_start <= pos_last_known_patho_ptc:
                comment = (
                    comment
                    + f" Frameshift variant starts upstream of p.His1184 and is predicted to lead to an alternative C-terminal end."
                )
                result = True
                strength = evidence_strength.STRONG
            elif transcript.var_start > pos_last_known_patho_ptc:
                comment = (
                    comment
                    + f" Frameshift variant starts downstream of p.Tyr1183 and is predicted to lead to an alternative C-terminal end."
                )
                result = True
                strength = evidence_strength.SUPPORTING
            else:
                raise ValueError(
                    f"Variant start in transcript {transcript.transcript_id} is {transcript.var_start} is neither smaller nor equal or bigger than 3552."
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

    @classmethod
    def assess_pvs1_splice_palb2(
        cls,
        transcript: TranscriptInfo_intronic,
        prediction_dict: dict[str, float],
        threshold: Threshold,
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        try:
            prediction_value = prediction_dict[threshold.name]
            prediction = assess_prediction_tool(threshold, prediction_value)
        except KeyError:
            prediction = None
        if not transcript.are_exons_skipped or prediction:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = (
                f"No splicing alteration predicted for {transcript.transcript_id}."
            )
        elif transcript.is_NMD and not transcript.is_reading_frame_preserved:
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = (
                f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            )
        elif not transcript.is_reading_frame_preserved:
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Reading frame is not preserved and truncated/altered region is disease relevant."
        else:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and reading frame is preserved."
            if transcript.is_truncated_region_disease_relevant:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = (
                    comment + f" Truncated region is critical to protein function."
                )
            else:
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                ):
                    result = True
                    strength = evidence_strength.STRONG
                    comment = (
                        comment
                        + f" Splicing alteration removes more than {threshold_diff_len_prot_percent} of coding sequence."
                    )
                elif (
                    abs(transcript.diff_len_protein_percent)
                    > threshold_diff_len_prot_percent
                ):
                    # In case of an increase in protein lenght, the diff_len_protein_percent is negative
                    result = True
                    strength = evidence_strength.MODERATE
                    comment = (
                        comment
                        + f" Splicing alteration inserts more than {threshold_diff_len_prot_percent} of coding sequence."
                    )
                else:
                    result = True
                    strength = evidence_strength.SUPPORTING
                    comment = (
                        comment
                        + f" Splicing alteration removes less then {threshold_diff_len_prot_percent} of coding sequence."
                    )
        return RuleResult(
            "PVS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
