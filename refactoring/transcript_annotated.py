#!/usr/bin/env python3

from __future__ import annotations

import pyensembl
from dataclasses import dataclass, field
from refactoring.genotoscope_protein_len_diff import calculate_prot_len_diff

from refactoring.variant import VariantInfo, TranscriptInfo
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
)
from refactoring.genotoscope_reading_frame_preservation import (
    assess_reading_frame_preservation,
)
from refactoring.genotoscope_exists_alternative_start_codon import (
    find_alternative_start_codon,
)


@dataclass
class TranscriptInfo_exonic(TranscriptInfo):
    """
    Class containing exonic variant specific annotation
    """

    ref_transcript: pyensembl.transcript.Transcript = field(init=True)
    var_seq: str = ""
    diff_len: int = 0
    diff_len_protein_percent: float = 0
    len_change_in_repetitive_region: bool = False
    truncated_exon_relevant: bool = False
    transcript_disease_relevant: bool = False
    is_NMD: bool = False
    is_reading_frame_preserved: bool = True

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
        is_NMD, NMD_affected_exons = assess_NMD_exonic_variant(
            transcript, variant, ref_transcript, var_seq, diff_len
        )
        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
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
            var_seq=var_seq,
            diff_len_protein_percent=diff_len_protein_percent,
            len_change_in_repetitive_region=False,
            is_NMD=is_NMD,
            is_reading_frame_preserved=is_reading_frame_preserved,
        )


@dataclass
class TranscriptInfo_intronic(TranscriptInfo):
    """
    Class containing intronic variant specific annotations
    """

    ref_transcript: pyensembl.transcript.Transcript
    diff_len: float
    diff_len_protein_percent: float
    are_exons_skipped: bool = False
    len_change_in_repetitive_region: bool = False
    is_NMD: bool = False
    transcript_disease_relevant: bool = False
    truncated_exon_relevant: bool = False
    is_reading_frame_preserved: bool = False

    @classmethod
    def annotate(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ) -> TranscriptInfo_intronic:
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
        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
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
            len_change_in_repetitive_region=False,
            is_NMD=is_NMD,
            transcript_disease_relevant=True,
            truncated_exon_relevant=True,
            is_reading_frame_preserved=is_reading_frame_preserved,
        )


@dataclass
class TranscriptInfo_start_loss(TranscriptInfo):
    """
    Class containing start loss specific annotation
    """

    ref_transcript: pyensembl.transcript.Transcript = field(init=True)
    exists_alternative_start_codon: bool
    position_alternative_start_codon: list[int]
    pathogenic_variant_between_start_and_stop: bool

    @classmethod
    def annotate(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ) -> TranscriptInfo_exonic:
        ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            transcript.transcript_id
        )
        diff_len, var_coding_seq = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        position_alternative_start_codon = find_alternative_start_codon(
            variant, ref_transcript, var_coding_seq
        )
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
            position_alternative_start_codon=position_alternative_start_codon,
            pathogenic_variant_between_start_and_stop=False,
        )
