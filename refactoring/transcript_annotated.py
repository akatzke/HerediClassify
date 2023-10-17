#!/usr/bin/env python3

from __future__ import annotations

import pyensembl
from dataclasses import dataclass, field
import pathlib

from refactoring.genotoscope_protein_len_diff import calculate_prot_len_diff

from refactoring.variant import VARTYPE_GROUPS, VariantInfo, TranscriptInfo
from refactoring.variant_annotate import ClinVar, Variant_annotated
from refactoring.genotoscope_exon_skipping import assess_exon_skipping
from refactoring.genotoscope_assess_NMD import (
    assess_NMD_exonic_variant,
    assess_NMD_intronic_variant,
)
from refactoring.genotoscope_construct_variant_sequence import (
    construct_variant_coding_seq_exonic_variant,
    construct_variant_coding_seq_intronic_variant,
)
from refactoring.genotoscope_assess_NMD import (
    assess_NMD_exonic_variant,
    assess_NMD_intronic_variant,
    assess_NMD_threshold
)
from refactoring.genotoscope_reading_frame_preservation import (
    assess_reading_frame_preservation,
)
from refactoring.genotoscope_exists_alternative_start_codon import (
    assess_alternative_start_codon,
)

from refactoring.genotoscope_protein_len_diff_repetitive_region import (
    check_prot_len_change_in_repetitive_region,
    check_prot_len_change_in_repetitive_region_start_loss,
    check_prot_len_change_in_repetitive_region_exon,
)

from refactoring.clinvar_region import (
    check_clinvar_NMD_exon,
    check_clinvar_start_alt_start,
    check_clinvar_truncated_region,
)



@dataclass
class TranscriptInfo_annot(TranscriptInfo):
    """
    Abstract class for all Transcript_Info annotation classes
    """

    ref_transcript: pyensembl.transcript.Transcript = field(init=True)
    diff_len_protein_percent: float = 0
    len_change_in_repetitive_region: bool = False
    is_truncated_region_disease_relevant: bool = False
    pathogenic_variants_truncated_region: list[str] = field(default_factory = list)


def annotate_transcripts(variant: Variant_annotated) -> list[TranscriptInfo_annot]:
    annotated_transcripts = []
    for transcript in variant.transcript_info:
        if transcript.var_type in VARTYPE_GROUPS.EXONIC.value:
            annotated_transcript = TranscriptInfo_exonic.annotate(variant, transcript)
        elif transcript.var_type in VARTYPE_GROUPS.INTRONIC.value:
            annotated_transcript = TranscriptInfo_intronic.annotate(variant, transcript)
        elif transcript.var_type in VARTYPE_GROUPS.START_LOST.value:
            annotated_transcript = TranscriptInfo_start_loss.annotate(variant, transcript)
        else:
            break
        annotated_transcripts.append(annotated_transcript)
    if len(annotated_transcripts) == 0:
        raise TypeError("No annotated transcripts created")
    return annotated_transcripts


def annotate_transcripts_acmg_specification(variant: Variant_annotated) -> list[TranscriptInfo_annot]:
    """
    Check if ACMG classification is sufficient to define complete variant interpretation
    """
    annotated_transcripts = []
    for transcript in variant.transcript_info:
        if transcript.var_type in VARTYPE_GROUPS.EXONIC.value:
            try:
                annotated_transcript = TranscriptInfo_exonic.annotate_acmg_specification(variant, transcript)
            except:
                annotated_transcript = TranscriptInfo_exonic.annotate(variant, transcript)
        elif transcript.var_type in VARTYPE_GROUPS.INTRONIC.value:
            try:
                annotated_transcript = TranscriptInfo_intronic.annotate_acmg_specification(variant, transcript)
            except:
                annotated_transcript = TranscriptInfo_intronic.annotate(variant, transcript)
        elif transcript.var_type in VARTYPE_GROUPS.START_LOST.value:
            try:
                annotated_transcript = TranscriptInfo_start_loss.annotate_acmg_specification(variant, transcript)
            except:
                annotated_transcript = TranscriptInfo_start_loss.annotate(variant, transcript)
        else:
            break
        annotated_transcripts.append(annotated_transcript)
    if len(annotated_transcripts) == 0:
        raise TypeError("No annotated transcripts created")
    return annotated_transcripts

@dataclass
class TranscriptInfo_exonic(TranscriptInfo_annot):
    """
    Class containing exonic variant specific annotation
    """

    is_NMD: bool = False
    is_reading_frame_preserved: bool = True

    @classmethod
    def get_arguments(cls) -> list:
        return []


    @classmethod
    def annotate(
            cls, variant: VariantInfo, transcript: TranscriptInfo
    ) -> TranscriptInfo_exonic:
        """
        Perform annotation for exonic variants
        """
        ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            transcript.transcript_id
        )
        var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        #if threshold_NMD:
        #    is_NMD, NMD_affected_exons = assess_NMD_threshold(transcript, threshold_NMD, clin_transcript)
        #else:
        is_NMD, NMD_affected_exons = assess_NMD_exonic_variant(transcript, variant, ref_transcript, var_seq, diff_len)
        #if functional_relevant_region:
        #    # Check if variant is located in defined functionally relevant region from ACMG guidelines
        #    is_truncated_exon_relevant, comment_truncated_exon_relevant = check_if_variant_affects_functional_region(variant, functional_relevant_region)
        #else:
        truncated_exon_ClinVar = check_clinvar_NMD_exon(variant, NMD_affected_exons)
        is_truncated_exon_relevant = truncated_exon_ClinVar.pathogenic
        comment_truncated_exon_relevant = truncated_exon_ClinVar.ids
        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
        if diff_len_protein_percent != 0:
            len_change_in_repetitive_region = check_prot_len_change_in_repetitive_region(variant, variant.genomic_start, variant.genomic_end, ref_transcript, pathlib.Path("/home/katzkean/variant_classification/genotoscope_data/uniprot/uniprot_feb_2021//home/katzkean/variant_classification/genotoscope_data/uniprot/uniprot_feb_2021/UP000005640_9606_repeat_hg19_16columns.bed"))
        else:
            len_change_in_repetitive_region = False
        return TranscriptInfo_exonic(
            transcript_id=transcript.transcript_id,
            var_type=transcript.var_type,
            var_hgvs=transcript.var_hgvs,
            var_start=transcript.var_start,
            var_stop=transcript.var_stop,
            var_protein=transcript.var_protein,
            exon=transcript.exon,
            intron=transcript.intron,
            ref_transcript=ref_transcript,
            diff_len_protein_percent=diff_len_protein_percent,
            len_change_in_repetitive_region=len_change_in_repetitive_region,
            is_NMD=is_NMD,
            is_reading_frame_preserved=is_reading_frame_preserved,
            is_truncated_region_disease_relevant=is_truncated_exon_relevant,
            pathogenic_variants_truncated_region = comment_truncated_exon_relevant
        )


@dataclass
class TranscriptInfo_intronic(TranscriptInfo_annot):
    """
    Class containing intronic variant specific annotations
    """

    are_exons_skipped: bool = False
    is_NMD: bool = False
    is_reading_frame_preserved: bool = False

    @classmethod
    def annotate(
            cls, variant: VariantInfo, transcript: TranscriptInfo
    ) -> TranscriptInfo_intronic:
        """
        Perform annotation specific for intronic variants
        """
        ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            transcript.transcript_id
        )
        (
            exons_skipped,
            are_exons_skipped,
            skipped_exon_start,
            skipped_exon_end,
            start_codon_exon_skipped,
            stop_codon_exon_skipped,
            coding_exon_skipped,
        ) = assess_exon_skipping(transcript, variant, ref_transcript)
        var_seq, diff_len = construct_variant_coding_seq_intronic_variant(
            transcript,
            variant,
            ref_transcript,
            skipped_exon_start,
            skipped_exon_end,
            are_exons_skipped,
            start_codon_exon_skipped,
            stop_codon_exon_skipped,
            coding_exon_skipped,
        )
        is_NMD, NMD_affected_exons = assess_NMD_intronic_variant(
            transcript,
            variant,
            ref_transcript,
            are_exons_skipped,
            exons_skipped,
            start_codon_exon_skipped,
            stop_codon_exon_skipped,
            coding_exon_skipped,
            var_seq,
            diff_len,
        )
        print(NMD_affected_exons)
        skipped_exon_ClinVar = check_clinvar_NMD_exon(variant, NMD_affected_exons)
        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
        if diff_len != 0:
            len_change_in_repetitive_region = check_prot_len_change_in_repetitive_region_exon(variant, ref_transcript, NMD_affected_exons, pathlib.Path("/home/katzkean/variant_classification/genotoscope_data/uniprot/uniprot_feb_2021/UP000005640_9606_repeat_hg19_16columns.bed"))
        else:
            len_change_in_repetitive_region = False
        return TranscriptInfo_intronic(
            transcript_id=transcript.transcript_id,
            var_type=transcript.var_type,
            var_hgvs=transcript.var_hgvs,
            var_start=transcript.var_start,
            var_stop=transcript.var_stop,
            var_protein=transcript.var_protein,
            exon=transcript.exon,
            intron=transcript.intron,
            ref_transcript=ref_transcript,
            diff_len_protein_percent=diff_len_protein_percent,
            are_exons_skipped=are_exons_skipped,
            len_change_in_repetitive_region=len_change_in_repetitive_region,
            is_NMD=is_NMD,
            is_truncated_region_disease_relevant=skipped_exon_ClinVar.pathogenic,
            pathogenic_variants_truncated_region=skipped_exon_ClinVar.ids,
            is_reading_frame_preserved=is_reading_frame_preserved,
        )


@dataclass
class TranscriptInfo_start_loss(TranscriptInfo_annot):
    """
    Class containing start loss specific annotation
    """

    exists_alternative_start_codon: bool = False
    position_alternative_start_codon: list[int] = field(default_factory=list)

    @classmethod
    def annotate(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ) -> TranscriptInfo_exonic:
        ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            transcript.transcript_id
        )
        var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        (
            exists_alternative_start_codon,
            position_alternative_start_codon,
        ) = assess_alternative_start_codon(variant, ref_transcript, var_seq)
        pathogenic_variants_between_start_and_alt_start = check_clinvar_start_alt_start(
            ref_transcript, variant, position_alternative_start_codon
        )
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
        if diff_len != 0:
            len_change_in_repetitive_region = check_prot_len_change_in_repetitive_region_start_loss(variant, ref_transcript, position_alternative_start_codon, pathlib.Path("/home/katzkean/variant_classification/genotoscope_data/uniprot/uniprot_feb_2021//home/katzkean/variant_classification/genotoscope_data/uniprot/uniprot_feb_2021/UP000005640_9606_repeat_hg19_16columns.bed"))
        else:
            len_change_in_repetitive_region = False
        return TranscriptInfo_start_loss(
            transcript_id=transcript.transcript_id,
            var_type=transcript.var_type,
            var_hgvs=transcript.var_hgvs,
            var_start=transcript.var_start,
            var_stop=transcript.var_stop,
            var_protein=transcript.var_protein,
            exon=transcript.exon,
            intron=transcript.intron,
            ref_transcript=ref_transcript,
            exists_alternative_start_codon=exists_alternative_start_codon,
            position_alternative_start_codon=position_alternative_start_codon,
            is_truncated_region_disease_relevant=pathogenic_variants_between_start_and_alt_start.pathogenic,
            pathogenic_variants_truncated_region=pathogenic_variants_between_start_and_alt_start.ids,
            diff_len_protein_percent=diff_len_protein_percent,
            len_change_in_repetitive_region=len_change_in_repetitive_region
        )
