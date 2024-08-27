#!/usr/bin/env python3

import pathlib

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from information import Classification_Info, Info
from variant import RNAData, TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from acmg_rules.functional_splicing_assay_utils import (
    adjust_strength_according_to_rna_data_pvs1,
)


class Pvs1(abstract_rule):
    """
    PVS1: Loss of function
    Devided into three separate parts: Frameshift, splice and start_loss
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
        threshold_diff_len_prot_percent: float,
        splice_assay: Optional[list[RNAData]],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        mane_path: pathlib.Path,
    ) -> RuleResult:
        results = {}
        for transcript in annotated_transcripts:
            if isinstance(transcript, TranscriptInfo_exonic):
                result_frameshift = Pvs1.assess_pvs1_frameshift_PTC(
                    transcript, threshold_diff_len_prot_percent
                )
                results[transcript.transcript_id] = result_frameshift
            elif isinstance(transcript, TranscriptInfo_intronic):
                result_splice = Pvs1.assess_pvs1_splice(
                    transcript,
                    prediction_dict,
                    threshold,
                    threshold_diff_len_prot_percent,
                )
                if splice_assay:
                    result_splice = adjust_strength_according_to_rna_data_pvs1(
                        splice_assay, result_splice
                    )
                results[transcript.transcript_id] = result_splice
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result_start_loss = Pvs1.assess_pvs1_start_loss(transcript)
                results[transcript.transcript_id] = result_start_loss
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
    def assess_pvs1_start_loss(
        cls, transcript: TranscriptInfo_start_loss
    ) -> RuleResult:
        """
        Assess PVS1 for start lost variants
        """
        if not transcript.exists_alternative_start_codon:
            comment = f"No alternative start codons were detected in transcript {transcript.transcript_id}."
            result = True
            strength = evidence_strength.MODERATE
        else:
            comment = f"Alternative start codon observed in transcript {transcript.transcript_id}."
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + f" Alternative start codon leads to the exclusion of a disease relevant protein region. "
                    + transcript.comment_truncated_region
                )
                result = True
                strength = evidence_strength.MODERATE
            else:
                comment = (
                    comment
                    + f" No pathogenic variant detected between start codon and alternative start codon."
                )
                result = True
                strength = evidence_strength.SUPPORTING
        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )

    @classmethod
    def assess_pvs1_start_loss_pathogenic_very_strong(cls) -> RuleResult:
        """
        Assess PVS1 for start loss variants, that are always PVS1 very strong
        """
        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            True,
            evidence_strength.VERY_STRONG,
            "For this gene, start loss variants always give PVS1 with very strong evidence strength.",
        )

    @classmethod
    def assess_pvs1_splice(
        cls,
        transcript: TranscriptInfo_intronic,
        prediction_dict: dict[str, float],
        threshold: Threshold,
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        """
        Assess PVS1 for splice variants
        """
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if not transcript.are_exons_skipped or not num_thresholds_met:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"No splicing alteration predicted for transcript {transcript.transcript_id}."
        elif transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} undergoes NMD."
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + f" Skipped exon is disease relevant. "
                    + transcript.comment_truncated_region
                )
                result = True
            else:
                comment = (
                    comment + " Skipped exon is not considered to be disease relevant."
                )
                result = False
            strength = evidence_strength.VERY_STRONG
        elif (
            transcript.are_exons_skipped
            and not transcript.is_NMD
            and not transcript.is_reading_frame_preserved
        ):
            comment = f"Transcript {transcript.transcript_id} does not undergo NMD and reading frame is not preserved."
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + f" Skipped exon is considered disease relevant. "
                    + transcript.comment_truncated_region
                )
                result = True
                strength = evidence_strength.STRONG
            else:
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                ):
                    comment = (
                        comment
                        + f" Protein length change of {transcript.diff_len_protein_percent} observed."
                    )
                    result = True
                    strength = evidence_strength.STRONG
                else:
                    comment = (
                        comment
                        + f" Protein length change of {transcript.diff_len_protein_percent} observed."
                    )
                    result = True
                    strength = evidence_strength.MODERATE
        elif transcript.are_exons_skipped and transcript.is_reading_frame_preserved:
            comment = f"Transcript {transcript.transcript_id} does not undergo NMD and reading frame is preserved."
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + f" Skipped exon is considered disease relevant region."
                    + transcript.comment_truncated_region
                )
                result = True
                strength = evidence_strength.STRONG
            else:
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                ):
                    comment = (
                        comment
                        + f" Protein length change of {transcript.diff_len_protein_percent} observed."
                    )
                    result = True
                    strength = evidence_strength.STRONG
                else:
                    comment = (
                        comment
                        + f" Protein length change of {transcript.diff_len_protein_percent} observed."
                    )
                    result = True
                    strength = evidence_strength.MODERATE
        else:
            comment = f"Transcript {transcript.transcript_id} does not fulfill any PVS1 splicing criteria."
            result = False
            strength = evidence_strength.VERY_STRONG
        return RuleResult(
            "PVS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )

    @classmethod
    def assess_pvs1_frameshift_PTC(
        cls, transcript: TranscriptInfo_exonic, threshold_diff_len_prot_percent: float
    ) -> RuleResult:
        """
        Assess PVS1 for frameshift variants
        """
        if transcript.is_NMD:
            comment = (
                f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            )
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + "Truncated region is disease relevant. "
                    + transcript.comment_truncated_region
                )
                result = True
            else:
                comment = comment + "Truncated region is not disease relevant."
                result = False
            strength = evidence_strength.VERY_STRONG
        else:
            comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD."
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + "Truncated region is disease relevant. "
                    + transcript.comment_truncated_region
                )
                result = True
                strength = evidence_strength.STRONG
            else:
                if (
                    transcript.diff_len_protein_percent
                    > threshold_diff_len_prot_percent
                ):
                    comment = (
                        comment
                        + f" Protein length change of {transcript.diff_len_protein_percent} observed."
                    )
                    result = True
                    strength = evidence_strength.STRONG
                else:
                    comment = (
                        comment
                        + f" Protein length change of {transcript.diff_len_protein_percent} observed."
                    )
                    result = True
                    strength = evidence_strength.MODERATE
        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
