#!/usr/bin/env python3

import pathlib

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
    assess_thresholds,
)
from acmg_rules.functional_splicing_assay_utils import (
    adjust_strength_according_to_rna_data_pvs1,
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
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
                class_info.SPLICE_RESULT,
                class_info.SPLICING_ASSAY,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        pos_last_known_patho_ptc_dict: dict[str, int],
        threshold_diff_len_prot_percent: float,
        splice_result: Optional[RuleResult],
        splice_assay: Optional[list[RNAData]],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        mane_path: pathlib.Path,
    ):
        results = {}
        for transcript in annotated_transcripts:
            if isinstance(transcript, TranscriptInfo_exonic):
                try:
                    pos_last_known_patho_ptc = pos_last_known_patho_ptc_dict[
                        transcript.transcript_id
                    ]
                except KeyError:
                    raise KeyError(
                        f"Transcript {transcript.transcript_id} not in disease relevant transcripts: {pos_last_known_patho_ptc_dict.keys()}. Transcript should have been filtered out earlier."
                    )
                result = cls.assess_pvs1_frameshift_PTC_brca1(
                    transcript,
                    pos_last_known_patho_ptc,
                )
                results[transcript.transcript_id] = result
            elif isinstance(transcript, TranscriptInfo_intronic):
                if splice_result is None:
                    splice_result = cls.assess_pvs1_splice_brca1(
                        transcript,
                        prediction_dict,
                        threshold,
                        threshold_diff_len_prot_percent,
                    )
                if splice_assay:
                    splice_result = adjust_strength_according_to_rna_data_pvs1(
                        splice_assay, splice_result
                    )
                results[transcript.transcript_id] = splice_result
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result = cls.assess_pvs1_start_loss_pathogenic_very_strong()
                results[transcript.transcript_id] = result
        if len(results) == 0:
            result = RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                comment=f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}.",
            )
        else:
            result = summarise_results_per_transcript(results, "PVS1", mane_path)
        return result

    @classmethod
    def assess_pvs1_frameshift_PTC_brca1(
        cls,
        transcript: TranscriptInfo_exonic,
        pos_last_known_patho_ptc: int,
    ) -> RuleResult:
        if transcript.is_NMD:
            comment = (
                f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            )
            if transcript.is_affected_exon_disease_relevant:
                comment = (
                    comment + " Variant is present in biologically relevant transcript."
                )
                result = True
            else:
                comment = (
                    comment
                    + " Variant is not present in biologically relevant transcript."
                )
                result = False
            strength = evidence_strength.VERY_STRONG
        else:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD."
            if transcript.ptc <= pos_last_known_patho_ptc:
                comment = (
                    comment
                    + " Truncated/altered region is critical to protein function."
                )
                result = True
                strength = evidence_strength.STRONG
            else:
                comment = (
                    comment
                    + " Role of truncated/alterend region in protein function is unknown."
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
    def assess_pvs1_splice_brca1(
        cls,
        transcript: TranscriptInfo_intronic,
        prediction_dict: dict[str, float],
        threshold: Threshold,
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if not transcript.are_exons_skipped and not num_thresholds_met:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"No splicing alteration predicted for transcript {transcript.transcript_id}."
        elif not transcript.coding_exon_skipped:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"Predicted alteration does not affect coding sequence in transcript {transcript.transcript_id}."
        elif transcript.start_codon_exon_skipped:
            result = True
            strength = evidence_strength.VERY_STRONG
            comment = f"Predicted alteration is non-coding (initation codon skipped) in transcript {transcript.transcript_id}."
        elif transcript.is_NMD:
            comment = (
                f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            )
            if transcript.is_affected_exon_disease_relevant:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = comment + (
                    " Skipped exon is present in biologically-relevant transcript."
                )
            else:
                result = False
                strength = evidence_strength.VERY_STRONG
                comment = (
                    comment
                    + " Skipped exon is absent in biologically-relevant transcript."
                )
        elif not transcript.is_reading_frame_preserved and not transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and reading frame is not preserved."
            if transcript.is_truncated_region_disease_relevant:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = (
                    comment
                    + f" Target region is critical to protein function (BRCT domain). "
                    + transcript.comment_truncated_region
                )
            else:
                comment = (
                    comment + " Exon is absent from biologically relevant transcript."
                )
                result = False
                strength = evidence_strength.VERY_STRONG
        elif transcript.is_reading_frame_preserved:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and reading frame is preserved."
            if transcript.is_truncated_region_disease_relevant:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = comment + f" Target region is critical to protein function."
            else:
                comment = comment + f" Role of target region is unkown."
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                ):
                    result = True
                    strength = evidence_strength.MODERATE
                    comment = (
                        comment
                        + f" Splicing alteration removes >{threshold_diff_len_prot_percent} of coding sequence."
                    )
                else:
                    result = True
                    strength = evidence_strength.SUPPORTING
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
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
