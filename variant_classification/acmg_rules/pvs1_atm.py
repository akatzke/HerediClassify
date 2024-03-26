#!/usr/bin/env python3

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from acmg_rules.pvs1 import Pvs1
from acmg_rules.pvs1_brca2 import Pvs1_brca2
from information import Classification_Info, Info
from variant import RNAData, TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from acmg_rules.functional_splicing_assay_utils import (
    adjust_strength_according_to_rna_data_pvs1,
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
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
                class_info.SPLICING_ASSAY,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        pos_last_known_patho_ptc_dict: dict[str, int],
        splice_result: Optional[RuleResult],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        splice_assay: Optional[list[RNAData]],
    ):
        results = []
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
                result = Pvs1_brca2.assess_pvs1_frameshift_PTC_brca2(
                    transcript, pos_last_known_patho_ptc
                )
                results.append(result)
            elif isinstance(transcript, TranscriptInfo_intronic):
                if splice_result is None:
                    splice_result = cls.assess_pvs1_splice_atm(
                        transcript, prediction_dict, threshold
                    )
                if splice_assay:
                    splice_result = adjust_strength_according_to_rna_data_pvs1(
                        splice_assay, splice_result
                    )
                results.append(splice_result)
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result = cls.assess_pvs1_start_loss_atm(transcript)
                results.append(result)
        if len(results) == 0:
            if not annotated_transcripts:
                comment = "No annotated transcripts provided, PVS1 can not be applied."
            else:
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
    def assess_pvs1_start_loss_atm(
        cls, transcript: TranscriptInfo_start_loss
    ) -> RuleResult:
        """
        Assess PVS1 for start lost variants
        """
        if transcript.is_truncated_region_disease_relevant:
            comment = (
                f"Alternative start codon leads to the exclusion of a disease relevant region. "
                + transcript.comment_truncated_region
            )
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

    @classmethod
    def assess_pvs1_splice_atm(
        cls,
        transcript: TranscriptInfo_intronic,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        disease_relevant_transcript = "ENST00000675843"
        if transcript.transcript_id != disease_relevant_transcript:
            raise ValueError(
                f"Transcript {transcript.transcript_id} is not disease relevant transcript {disease_relevant_transcript}. This is hard coded and not changable. Transcript should have been filterd out earlier."
            )
        if not transcript.are_exons_skipped or not num_thresholds_met:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = (
                f"No splicing alteration predicted for {transcript.transcript_id}."
            )
        elif not transcript.coding_exon_skipped:
            result = False
            strength = evidence_strength.VERY_STRONG
            comment = f"Exon skipping or use of cryptic splice site does not affect the coding sequence in transcript {transcript.transcript_id}."
        elif (
            transcript.affected_exon["exon_no"] >= 2
            and transcript.affected_exon["exon_no"] <= 38
        ):
            if not transcript.is_reading_frame_preserved:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            else:
                result = True
                strength = evidence_strength.STRONG
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and ."
        elif (
            transcript.affected_exon["exon_no"] >= 39
            and transcript.affected_exon["exon_no"] <= 63
        ):
            if not transcript.is_reading_frame_preserved and transcript.is_NMD:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = f"Transcript {transcript.transcript_id} is predicted to undergo NMD."
            elif not transcript.is_reading_frame_preserved and not transcript.is_NMD:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and reading frame is not preserved. Truncated region is disease relevant."
            elif transcript.is_reading_frame_preserved:
                result = True
                strength = evidence_strength.VERY_STRONG
                comment = f"Transcript {transcript.transcript_id} is not predicted to undergo NMD and reading frame is preserved. Truncated region is disease relevant."
            else:
                result = False
                strength = evidence_strength.VERY_STRONG
                comment = f"Splicing alteration in transcript {transcript.transcript_id} not predicted to be pathogenic."
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
