#!/usr/bin/env python3

from __future__ import annotations

import pyensembl
import pathlib
import logging
from dataclasses import dataclass, field
from collections.abc import Callable

from variant import VariantInfo, TranscriptInfo, Variant
from ensembl import ensembl
from genotoscope_construct_variant_sequence import (
    construct_variant_coding_seq_exonic_variant,
    construct_variant_coding_seq_intronic_variant,
)
from genotoscope_assess_NMD import (
    assess_NMD_exonic_variant,
    assess_NMD_intronic_variant,
)
from genotoscope_reading_frame_preservation import (
    assess_reading_frame_preservation,
)
from genotoscope_protein_len_diff import (
    calculate_prot_len_diff,
    calculate_prot_len_diff_start_loss,
)
from variant_in_critical_region import (
    check_variant_in_critical_region_exon,
)
from clinvar_region import (
    check_clinvar_NMD_exon,
    check_clinvar_start_alt_start,
    check_clinvar_region,
)
from genotoscope_exon_skipping import assess_exon_skipping
from genotoscope_protein_len_diff_repetitive_region import (
    check_prot_len_change_in_repetitive_region_exon,
)
from genotoscope_exists_alternative_start_codon import (
    assess_alternative_start_codon,
)
from utils import (
    check_intersection_with_bed,
    check_bed_intersect_start_loss,
)
from var_type import VARTYPE_GROUPS
from information import Classification_Info, Info


logger = logging.getLogger("GenOtoScope_Classify.config_annotation")


@dataclass
class TranscriptInfo_annot(TranscriptInfo):
    """
    Abstract class for all Transcript_Info annotation classes
    """

    ref_transcript: pyensembl.transcript.Transcript = field(init=True)
    diff_len_protein_percent: float = 0
    len_change_in_repetitive_region: bool = False
    is_truncated_region_disease_relevant: bool = False
    comment_truncated_region: str = ""


def annotate_transcripts(
    variant: Variant,
    fun_dict: dict[VARTYPE_GROUPS, Callable[[TranscriptInfo], TranscriptInfo_annot]],
) -> list[TranscriptInfo_annot]:
    annotated_transcripts = []
    for transcript in variant.transcript_info:
        if any(
            var_type in VARTYPE_GROUPS.EXONIC.value for var_type in transcript.var_type
        ):
            annot_fun = fun_dict[VARTYPE_GROUPS.EXONIC]
        elif any(
            var_type in VARTYPE_GROUPS.INTRONIC.value
            for var_type in transcript.var_type
        ):
            annot_fun = fun_dict[VARTYPE_GROUPS.INTRONIC]
        elif any(
            var_type in VARTYPE_GROUPS.START_LOST.value
            for var_type in transcript.var_type
        ):
            annot_fun = fun_dict[VARTYPE_GROUPS.INTRONIC]
        else:
            break
        try:
            annotated_transcript = annot_fun(transcript)
            annotated_transcripts.append(annotated_transcript)
        except Exception:
            logger.warning(f"Annotation function {annot_fun.func} did not execute.")
    return annotated_transcripts


def annotate_transcripts_acmg_specification(
    variant: Variant,
    fun_dict: dict[VARTYPE_GROUPS, dict[str, Callable[[], TranscriptInfo_annot]]],
) -> list[TranscriptInfo_annot]:
    """
    Check if ACMG classification is sufficient to define complete variant interpretation
    """
    annotated_transcripts = []
    for transcript in variant.transcript_info:
        if transcript.var_type in VARTYPE_GROUPS.EXONIC.value:
            try:
                annot_fun = fun_dict[VARTYPE_GROUPS.EXONIC]["acmg"]
                annotated_transcript = annot_fun()
            except Exception:
                annot_fun = fun_dict[VARTYPE_GROUPS.EXONIC]["general"]
                annotated_transcript = annot_fun()
        elif transcript.var_type in VARTYPE_GROUPS.INTRONIC.value:
            try:
                annot_fun = fun_dict[VARTYPE_GROUPS.INTRONIC]["acmg"]
                annotated_transcript = annot_fun()
            except Exception:
                annot_fun = fun_dict[VARTYPE_GROUPS.INTRONIC]["general"]
                annotated_transcript = annot_fun()
        elif transcript.var_type in VARTYPE_GROUPS.START_LOST.value:
            try:
                annot_fun = fun_dict[VARTYPE_GROUPS.START_LOST]["acmg"]
                annotated_transcript = annot_fun()
            except Exception:
                annot_fun = fun_dict[VARTYPE_GROUPS.START_LOST]["general"]
                annotated_transcript = annot_fun()
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
    def get_annotate(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.annotate,
            (
                class_info.VARIANT,
                class_info.CLINVAR_PATH,
                class_info.UNIPROT_REP_REGION_PATH,
                class_info.CRITICAL_REGION_PATH,
            ),
        )

    @classmethod
    def annotate(
        cls,
        variant: VariantInfo,
        path_clinvar: pathlib.Path,
        path_uniprot_rep: pathlib.Path,
        path_critical_region: pathlib.Path,
        transcript: TranscriptInfo,
    ) -> TranscriptInfo_exonic:
        """
        Perform annotation for exonic variants
        """
        ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
        var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        # if threshold_NMD:
        #    is_NMD, NMD_affected_exons = assess_NMD_threshold(transcript, threshold_NMD, clin_transcript)
        # else:
        is_NMD, NMD_affected_exons = assess_NMD_exonic_variant(
            transcript, variant, ref_transcript, var_seq, diff_len
        )
        if path_critical_region is not None:
            # Check if variant is located in defined functionally relevant region from ACMG guidelines
            if is_NMD:
                is_truncated_exon_relevant = check_variant_in_critical_region_exon(
                    variant, ref_transcript, NMD_affected_exons, path_critical_region
                )
                comment_truncated_exon_relevant = ""
            else:
                is_truncated_exon_relevant = check_intersection_with_bed(
                    variant,
                    variant.genomic_start,
                    variant.genomic_end,
                    ref_transcript,
                    path_critical_region,
                )
                comment_truncated_exon_relevant = ""
        else:
            if is_NMD:
                truncated_exon_ClinVar = check_clinvar_NMD_exon(
                    variant, NMD_affected_exons, path_clinvar
                )
                is_truncated_exon_relevant = truncated_exon_ClinVar.pathogenic
                comment_truncated_exon_relevant = f"The following relevant ClinVar are (likely) pathogenic: {truncated_exon_ClinVar.ids}"
            else:
                truncated_exon_ClinVar = check_clinvar_region(
                    variant, variant.genomic_start, variant.genomic_end, path_clinvar
                )
                is_truncated_exon_relevant = truncated_exon_ClinVar.pathogenic
                comment_truncated_exon_relevant = f"The following relevant ClinVar are (likely) pathogenic: {truncated_exon_ClinVar.ids}"

        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
        if diff_len_protein_percent != 0:
            len_change_in_repetitive_region = check_intersection_with_bed(
                variant,
                variant.genomic_start,
                variant.genomic_end,
                ref_transcript,
                path_uniprot_rep,
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
            diff_len_protein_percent=diff_len_protein_percent,
            len_change_in_repetitive_region=len_change_in_repetitive_region,
            is_NMD=is_NMD,
            is_reading_frame_preserved=is_reading_frame_preserved,
            is_truncated_region_disease_relevant=is_truncated_exon_relevant,
            comment_truncated_region=comment_truncated_exon_relevant,
        )

    @classmethod
    def get_annotate_acmg(cls, class_info: Classification_Info):
        class_info = class_info
        pass

    @classmethod
    def annotate_acmg_specification(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ):
        """
        This function implements gene specifications from ACMG for frameshift and other variants
        - NMD threshold
        - Functionally important regions
        """
        variant = variant
        transcript = transcript
        pass


@dataclass
class TranscriptInfo_intronic(TranscriptInfo_annot):
    """
    Class containing intronic variant specific annotations
    """

    are_exons_skipped: bool = False
    is_NMD: bool = False
    is_reading_frame_preserved: bool = False

    @classmethod
    def get_annotate(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.annotate,
            (
                class_info.VARIANT,
                class_info.CLINVAR_PATH,
                class_info.UNIPROT_REP_REGION_PATH,
                class_info.CRITICAL_REGION_PATH,
            ),
        )

    @classmethod
    def annotate(
        cls,
        variant: VariantInfo,
        path_clinvar: pathlib.Path,
        path_uniprot_rep: pathlib.Path,
        path_critical_region,
        transcript: TranscriptInfo,
    ) -> TranscriptInfo_intronic:
        """
        Perform annotation specific for intronic variants
        """
        ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
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
        if path_critical_region is not None:
            is_truncated_region_disease_relevant = (
                check_variant_in_critical_region_exon(
                    variant,
                    ref_transcript,
                    NMD_affected_exons,
                    path_critical_region,
                )
            )
            comment_truncated_region_disease_relevant = ""
        else:
            skipped_exon_ClinVar = check_clinvar_NMD_exon(
                variant, NMD_affected_exons, path_clinvar
            )
            is_truncated_region_disease_relevant = skipped_exon_ClinVar.pathogenic
            comment_truncated_region_disease_relevant = f"The following relevant ClinVar are (likely) pathogenic: {skipped_exon_ClinVar.ids}"
        is_reading_frame_preserved = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent = calculate_prot_len_diff(ref_transcript, var_seq)
        if diff_len != 0:
            len_change_in_repetitive_region = (
                check_prot_len_change_in_repetitive_region_exon(
                    variant, ref_transcript, NMD_affected_exons, path_uniprot_rep
                )
            )
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
            is_truncated_region_disease_relevant=is_truncated_region_disease_relevant,
            comment_truncated_region=comment_truncated_region_disease_relevant,
            is_reading_frame_preserved=is_reading_frame_preserved,
        )

    @classmethod
    def get_annotate_acmg(cls):
        pass

    @classmethod
    def annotate_acmg_specification(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ):
        """
        This function implements gene specifications from ACMG for frameshift and other variants
        - Preimplemented splice decision trees
        """
        variant = variant
        transcript = transcript
        pass


@dataclass
class TranscriptInfo_start_loss(TranscriptInfo_annot):
    """
    Class containing start loss specific annotation
    """

    exists_alternative_start_codon: bool = False
    position_alternative_start_codon: list[int] = field(default_factory=list)

    @classmethod
    def get_annotate(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.annotate,
            (
                class_info.VARIANT,
                class_info.CLINVAR_PATH,
                class_info.UNIPROT_REP_REGION_PATH,
                class_info.CRITICAL_REGION_PATH,
            ),
        )

    @classmethod
    def annotate(
        cls,
        variant: VariantInfo,
        path_clinvar: pathlib.Path,
        path_uniprot_rep: pathlib.Path,
        path_critical_region: pathlib.Path,
        transcript: TranscriptInfo,
    ) -> TranscriptInfo_start_loss:
        ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
        var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        (
            exists_alternative_start_codon,
            position_alternative_start_codon_genomic,
            position_alternative_start_codon_cDNA,
        ) = assess_alternative_start_codon(variant, ref_transcript, var_seq)
        if path_critical_region is not None:
            is_truncated_region_disease_relevant = check_bed_intersect_start_loss(
                variant,
                ref_transcript,
                position_alternative_start_codon_genomic,
                path_critical_region,
            )
            comment_truncated_region = ""
        else:
            pathogenic_variants_between_start_and_alt_start = (
                check_clinvar_start_alt_start(
                    ref_transcript,
                    variant,
                    position_alternative_start_codon_genomic,
                    path_clinvar,
                )
            )
            is_truncated_region_disease_relevant = (
                pathogenic_variants_between_start_and_alt_start.pathogenic
            )
            comment_truncated_region = f"The following relevant ClinVar entries are (likely) pathogenic: {pathogenic_variants_between_start_and_alt_start.ids}"
        diff_len_protein_percent = calculate_prot_len_diff_start_loss(
            ref_transcript, position_alternative_start_codon_cDNA
        )
        if diff_len != 0:
            len_change_in_repetitive_region = check_bed_intersect_start_loss(
                variant,
                ref_transcript,
                position_alternative_start_codon_genomic,
                path_uniprot_rep,
            )
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
            position_alternative_start_codon=position_alternative_start_codon_genomic,
            is_truncated_region_disease_relevant=is_truncated_region_disease_relevant,
            comment_truncated_region=comment_truncated_region,
            diff_len_protein_percent=diff_len_protein_percent,
            len_change_in_repetitive_region=len_change_in_repetitive_region,
        )

    @classmethod
    def get_annotate_acmg(cls, class_info: Classification_Info):
        class_info = class_info
        pass

    @classmethod
    def annotate_acmg_specification(
        cls, variant: VariantInfo, transcript: TranscriptInfo
    ):
        """
        This function implements gene specifications from ACMG for start loss variants
        """
        variant = variant
        transcript = transcript
        pass
