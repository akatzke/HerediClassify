#!/usr/bin/env python3

from typing import Callable

from variant_classification.acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
    summarise_results_per_transcript,
)
from variant_classification.information import Classification_Info, Info
from variant_classification.variant import TranscriptInfo, VariantInfo
from variant_classification.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
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
            (class_info.ANNOTATED_TRANSCRIPT_LIST, class_info.VARIANT),
        )

    @classmethod
    def assess_rule(
        cls, annotated_transcripts: list[TranscriptInfo], variant: VariantInfo
    ) -> RuleResult:
        results = []
        for transcript in annotated_transcripts:
            if type(transcript) is TranscriptInfo_exonic:
                result_frameshift = cls.assess_pvs1_frameshift_PTC(transcript)
                results.append(result_frameshift)
            elif type(transcript) is TranscriptInfo_intronic:
                result_splice = cls.assess_pvs1_splice(transcript)
                results.append(result_splice)
            elif type(transcript) is TranscriptInfo_start_loss:
                result_start_loss = cls.assess_pvs1_start_loss(transcript)
                results.append(result_start_loss)
        if len(results) == 0:
            comment = f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {variant.var_type}."
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
    def assess_pvs1_start_loss(
        cls, transcript: TranscriptInfo_start_loss
    ) -> RuleResult:
        """
        Assess PVS1 for start lost variants
        """
        if not transcript.exists_alternative_start_codon:
            comment = f"An alternative start code in transcript {transcript.transcript_id} at {transcript.position_alternative_start_codon} detected."
            result = False
            strength = evidence_strength.VERY_STRONG
        else:
            comment = f"No alternative start codon detected for transcript {transcript.transcript_id}."
            if transcript.is_truncated_region_disease_relevant:
                comment = f"Pathogenic variants detected between start codon and alternative start codon detected. \n ClinVar ID: {transcript.pathogenic_variants_truncated_region}"
                result = True
                strength = evidence_strength.MODERATE
            else:
                comment = "No pathogenic variant detected between start codon and alternative start codon."
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
    def assess_pvs1_splice(cls, transcript: TranscriptInfo_intronic) -> RuleResult:
        """
        Assess PVS1 for splice variants
        """
        if transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} undergoes NMD."
            if transcript.is_truncated_region_disease_relevant:
                comment = f"Skipped exon contais (likely) pathogenic variants and can therefore be considexoneredexon to be disease relevant. \n ClinVar ID: {transcript.pathogenic_variants_truncated_region}"
                result = True
            else:
                comment = "Skipped exon contains no (likely) pathogenic variants and is therefore not considered disease relevant."
                result = False
            strength = evidence_strength.VERY_STRONG
        elif (
            transcript.are_exons_skipped
            and not transcript.is_NMD
            and not transcript.is_reading_frame_preserved
        ):
            comment = f"Transcript {transcript.transcript_id} does not undergo NMD."
            if transcript.is_truncated_region_disease_relevant:
                comment = f"Skipped exon contais (likely) pathogenic variants and can therefore be considered to be disease relevant. \n ClinVar ID: {transcript.pathogenic_variants_truncated_region}"
                result = True
                strength = evidence_strength.STRONG
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Protein length change of"
                    result = True
                    strength = evidence_strength.STRONG
                else:
                    comment = "Something"
                    result = True
                    strength = evidence_strength.MODERATE
        elif transcript.are_exons_skipped and transcript.is_reading_frame_preserved:
            if transcript.is_truncated_region_disease_relevant:
                comment = "Something"
                result = True
                strength = evidence_strength.STRONG
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    result = True
                    strength = evidence_strength.STRONG
                else:
                    comment = "Something"
                    result = True
                    strength = evidence_strength.MODERATE
        else:
            comment = "Something"
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
        cls, transcript: TranscriptInfo_exonic
    ) -> RuleResult:
        """
        Assess PVS1 for frameshift variants
        """
        if transcript.is_NMD:
            if transcript.is_truncated_region_disease_relevant:
                comment = "Something"
                result = True
            else:
                comment = "Something"
                result = False
            strength = evidence_strength.VERY_STRONG
        else:
            if transcript.is_truncated_region_disease_relevant:
                comment = "Something"
                result = True
                strength = evidence_strength.STRONG
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    result = True
                    strength = evidence_strength.STRONG
                else:
                    comment = "Something"
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
