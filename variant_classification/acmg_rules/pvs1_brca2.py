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
from variant import TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)


class Pvs1_brca2(Pvs1):
    """
    PVS1: loss of function
    Following VEP guidelines for BRCA2
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
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
                class_info.POS_LAST_KNOWN_PATHO_PTC,
                class_info.SPLICE_RESULT,
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
                result = cls.assess_pvs1_frameshift_PTC_brca2(
                    transcript, pos_last_known_patho_ptc
                )
                results.append(result)
            elif isinstance(transcript, TranscriptInfo_intronic):
                if splice_result is None:
                    result = cls.assess_pvs1_splice_brca2(
                        transcript, threshold_diff_len_prot_percent
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
    def assess_pvs1_frameshift_PTC_brca2(
        cls, transcript: TranscriptInfo_exonic, pos_last_known_patho_ptc: int
    ) -> RuleResult:
        if transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD and in a disease relevant transcript."
            result = True
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

    @classmethod
    def assess_pvs1_splice_brca2(
        cls, transcript: TranscriptInfo_intronic, threshold_diff_len_prot_percent: float
    ) -> RuleResult:
        if not transcript.are_exons_skipped:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"No splicing alteration predicted for transcript {transcript.transcript_id}."
        if not transcript.coding_exon_skipped:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"Predicted alteration does not affect coding sequence {transcript.transcript_id}."
        elif transcript.start_codon_exon_skipped:
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Predicted alteration is non-coding (initation codon skipped) {transcript.transcript_id}."
        elif transcript.is_NMD:
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD. All exons in transcript are disease relevant."
        elif not transcript.is_reading_frame_preserved and not transcript.is_NMD:
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD. Reading frame is not preserved and in BRCA2 all regions are considered disease relevant."
        elif transcript.is_reading_frame_preserved:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and reading frame is not preserved."
            if transcript.is_truncated_region_disease_relevant:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = (
                    comment
                    + f" Transcript {transcript.transcript_id} is not predicted to undergo NMD. Reading frame is not preserved and in BRCA2 all regions are considered disease relevant."
                )
            else:
                comment = comment + " Role of target region is unknown."
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                ):
                    comment = (
                        comment
                        + f" Splicing alteration removes >{threshold_diff_len_prot_percent} of coding sequence."
                    )
                    if transcript.exon_skipped == 11:
                        result = True
                        strength = evidence_strength.VERY_STRONG
                        comment = comment + f" Exon 11 is skipped."
                    else:
                        result = False
                        strength = evidence_strength.VERY_STRONG
                        comment = comment + f" Exon 11 is not skipped."
                else:
                    result = False
                    strength = evidence_strength.VERY_STRONG
                    comment = (
                        comment
                        + f" Splicing alteration removes <{threshold_diff_len_prot_percent} of coding sequence."
                    )
        else:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"Splicing alteration in transcript {transcript.transcript_id} not predicted to be pathogenic."
        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
