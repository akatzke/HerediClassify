#!/usr/bin/env python3

from __future__ import annotations

import pyensembl
from dataclasses import dataclass, field
from refactoring.genotoscope_protein_len_diff import calculate_prot_len_diff

from refactoring.variant import VariantInfo, TranscriptInfo
from refactoring.variant_annotate import ClinVar
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
    assess_alternative_start_codon,
)
from refactoring.clinvar_region import (
    check_clinvar_NMD_exon,
    check_clinvar_start_alt_start,
    check_clinvar_truncated_region,
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
    is_NMD: bool = False
    is_truncated_exon_relevant: bool = False
    pathogenic_variants_truncated_exons: list[str] = field(default_factory=list)
    is_NMD_exon_relevant: bool = False
    pathogenic_variants_NMD_exon: list[str] = field(default_factory=list)
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
        print(NMD_affected_exons)
        if is_NMD:
            NMD_exon_ClinVar = check_clinvar_NMD_exon(variant, NMD_affected_exons)
            truncated_region_ClinVar = ClinVar(
                pathogenic=False, type="region", highest_classification=None
            )
        else:
            NMD_exon_ClinVar = ClinVar(
                pathogenic=False, type="region", highest_classification=None
            )
            truncated_region_ClinVar = check_clinvar_truncated_region(
                variant, ref_transcript
            )
        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
        if diff_len_protein_percent != 0:
            len_change_in_repetitive_region = (
                check_prot_len_change_in_repetitive_region(variant)
            )
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
            var_seq=var_seq,
            diff_len_protein_percent=diff_len_protein_percent,
            len_change_in_repetitive_region=len_change_in_repetitive_region,
            is_NMD=is_NMD,
            is_reading_frame_preserved=is_reading_frame_preserved,
            is_truncated_exon_relevant=truncated_region_ClinVar.pathogenic,
            pathogenic_variants_truncated_exons=truncated_region_ClinVar.ids,
            is_NMD_exon_relevant=NMD_exon_ClinVar.pathogenic,
            pathogenic_variants_NMD_exon=NMD_exon_ClinVar.ids,
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
    skipped_exon_relevant: bool = False
    pathogenic_variants_skipped_exon: list[str] = field(default_factory=list)
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
        print(NMD_affected_exons)
        skipped_exon_ClinVar = check_clinvar_NMD_exon(variant, NMD_affected_exons)
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
            diff_len=diff_len,
            skipped_exon_relevant=skipped_exon_ClinVar.pathogenic,
            pathogenic_variants_skipped_exon=skipped_exon_ClinVar.ids,
            is_reading_frame_preserved=is_reading_frame_preserved,
        )


@dataclass
class TranscriptInfo_start_loss(TranscriptInfo):
    """
    Class containing start loss specific annotation
    """

    exists_alternative_start_codon: bool = False
    position_alternative_start_codon: list[int] = field(default_factory=list)
    are_pathogenic_variants_between_start_and_alt_start: bool = False
    pathogenic_variants_between_start_and_alt_start: list[str] = field(
        default_factory=list
    )

    @classmethod
    def annotate(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ) -> TranscriptInfo_exonic:
        ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            transcript.transcript_id
        )
        var_coding_seq, diff_len = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        (
            exists_alternative_start_codon,
            position_alternative_start_codon,
        ) = assess_alternative_start_codon(variant, ref_transcript, var_coding_seq)
        pathogenic_variants_between_start_and_alt_start = check_clinvar_start_alt_start(
            ref_transcript, variant, position_alternative_start_codon
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
            exists_alternative_start_codon=exists_alternative_start_codon,
            position_alternative_start_codon=position_alternative_start_codon,
            are_pathogenic_variants_between_start_and_alt_start=pathogenic_variants_between_start_and_alt_start.pathogenic,
            pathogenic_variants_between_start_and_alt_start=pathogenic_variants_between_start_and_alt_start.ids,
        )
