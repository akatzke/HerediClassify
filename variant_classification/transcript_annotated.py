#!/usr/bin/env python3

from __future__ import annotations

import pyensembl
import pathlib
import logging
from dataclasses import dataclass, field
from collections.abc import Callable
from typing import Optional

from variant import VariantInfo, TranscriptInfo, Variant
from ensembl import ensembl
from genotoscope_construct_variant_sequence import (
    construct_variant_coding_seq_exonic_variant,
    construct_variant_coding_seq_intronic_variant,
)
from genotoscope_assess_NMD import (
    assess_NMD_threshold,
    assess_NMD_exonic_variant,
    assess_NMD_intronic_variant,
    get_affected_exon,
)
from genotoscope_reading_frame_preservation import (
    assess_reading_frame_preservation,
)
from genotoscope_protein_len_diff import (
    calculate_prot_len_diff,
    calculate_prot_len_diff_start_loss,
)
from custom_exceptions import (
    Pyensembl_no_coding_sequence,
    Pyensembl_transcript_not_found,
    Not_disease_relevant_transcript,
)
from variant_in_critical_region import (
    check_variant_in_critical_region_exon,
)
from clinvar_region import (
    check_clinvar_NMD_exon,
    check_clinvar_start_alt_start,
    check_clinvar_truncated_region,
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
from check_exon_disease_relevant import check_exon_disease_relevant
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
            annot_fun = fun_dict[VARTYPE_GROUPS.START_LOST]
        else:
            continue
        try:
            annotated_transcript = annot_fun(transcript)
            annotated_transcripts.append(annotated_transcript)
        except Pyensembl_no_coding_sequence:
            logger.warning(
                f"No coding sequence in pyensembl for transcript {transcript.transcript_id}.",
                exc_info=True,
            )
        except Pyensembl_transcript_not_found:
            logger.warning(
                f"No transcript in pyensembl for transcript {transcript.transcript_id}.",
                exc_info=True,
            )
    return annotated_transcripts


@dataclass
class TranscriptInfo_exonic(TranscriptInfo_annot):
    """
    Class containing exonic variant specific annotation
    """

    is_NMD: bool = False
    affected_exon: dict = field(default_factory=dict)
    is_reading_frame_preserved: bool = True
    frameshift: int = 0
    ptc: int = -1
    is_affected_exon_disease_relevant: bool = True

    @classmethod
    def get_annotate(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.annotate,
            (
                class_info.VARIANT,
                class_info.CLINVAR_PATH,
                class_info.CLINVAR_PATH_INDEL,
                class_info.UNIPROT_REP_REGION_PATH,
                class_info.CRITICAL_REGION_PATH,
                class_info.DISEASE_IRRELEVANT_EXONS_PATH,
                class_info.THRESHOLD_NMD,
            ),
        )

    @classmethod
    def annotate(
        cls,
        variant: VariantInfo,
        path_clinvar_snv: pathlib.Path,
        path_clinvar_indel: pathlib.Path,
        path_uniprot_rep: pathlib.Path,
        path_critical_region: Optional[pathlib.Path],
        path_disease_irrelevant_exons: Optional[pathlib.Path],
        nmd_threshold_dict: Optional[dict[str, int]],
        transcript: TranscriptInfo,
    ) -> TranscriptInfo_exonic:
        """
        Perform annotation for exonic variants
        """
        try:
            ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
        except ValueError:
            raise Pyensembl_transcript_not_found
        try:
            ref_transcript.coding_sequence
        except ValueError or AttributeError:
            raise Pyensembl_no_coding_sequence
        var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
            transcript, variant, ref_transcript
        )
        if nmd_threshold_dict is not None:
            try:
                nmd_threshold = nmd_threshold_dict[transcript.transcript_id]
            except KeyError:
                raise Not_disease_relevant_transcript
            is_NMD, NMD_affected_exons = assess_NMD_threshold(
                transcript, variant, ref_transcript, diff_len, nmd_threshold
            )
        else:
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
                    variant, NMD_affected_exons, path_clinvar_snv, path_clinvar_indel
                )
                is_truncated_exon_relevant = truncated_exon_ClinVar.pathogenic
                comment_truncated_exon_relevant = f"The following relevant ClinVar are (likely) pathogenic: {truncated_exon_ClinVar.ids}"
            else:
                truncated_exon_ClinVar = check_clinvar_truncated_region(
                    variant,
                    ref_transcript,
                    path_clinvar_snv,
                    path_clinvar_indel,
                )
                is_truncated_exon_relevant = truncated_exon_ClinVar.pathogenic
                comment_truncated_exon_relevant = f"The following relevant ClinVar are (likely) pathogenic: {truncated_exon_ClinVar.ids}"
        if not NMD_affected_exons:
            affected_exon = get_affected_exon(
                ref_transcript, transcript, variant, diff_len
            )[0]
        else:
            affected_exon = NMD_affected_exons[0]
        if is_NMD and path_disease_irrelevant_exons is not None:
            is_affected_exon_disease_relevant = check_exon_disease_relevant(
                path_disease_irrelevant_exons, NMD_affected_exons
            )
        else:
            is_affected_exon_disease_relevant = True
        is_reading_frame_preserved, frameshift = assess_reading_frame_preservation(
            diff_len
        )
        diff_len_protein_percent, ptc = calculate_prot_len_diff(
            ref_transcript, var_seq, diff_len
        )
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
            affected_exon=affected_exon,
            is_reading_frame_preserved=is_reading_frame_preserved,
            frameshift=frameshift,
            is_truncated_region_disease_relevant=is_truncated_exon_relevant,
            is_affected_exon_disease_relevant=is_affected_exon_disease_relevant,
            comment_truncated_region=comment_truncated_exon_relevant,
            ptc=ptc,
        )


@dataclass
class TranscriptInfo_intronic(TranscriptInfo_annot):
    """
    Class containing intronic variant specific annotations
    """

    are_exons_skipped: bool = False
    coding_exon_skipped: bool = False
    start_codon_exon_skipped: bool = False
    is_NMD: bool = False
    ptc: int = -1
    affected_exon: dict = field(default_factory=dict)
    is_reading_frame_preserved: bool = False
    is_affected_exon_disease_relevant: bool = True

    @classmethod
    def get_annotate(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.annotate,
            (
                class_info.VARIANT,
                class_info.CLINVAR_PATH,
                class_info.CLINVAR_PATH_INDEL,
                class_info.UNIPROT_REP_REGION_PATH,
                class_info.DISEASE_IRRELEVANT_EXONS_PATH,
                class_info.CRITICAL_REGION_PATH,
            ),
        )

    @classmethod
    def annotate(
        cls,
        variant: VariantInfo,
        path_clinvar_snv: pathlib.Path,
        path_clinvar_indel: pathlib.Path,
        path_uniprot_rep: pathlib.Path,
        path_disease_irrelevant_exons: Optional[pathlib.Path],
        path_critical_region: Optional[pathlib.Path],
        transcript: TranscriptInfo,
    ) -> TranscriptInfo_intronic:
        """
        Perform annotation specific for intronic variants
        """
        try:
            ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
        except ValueError:
            raise Pyensembl_transcript_not_found
        try:
            ref_transcript.coding_sequence
        except ValueError or AttributeError:
            raise Pyensembl_no_coding_sequence
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
                variant, NMD_affected_exons, path_clinvar_snv, path_clinvar_indel
            )
            is_truncated_region_disease_relevant = skipped_exon_ClinVar.pathogenic
            comment_truncated_region_disease_relevant = f"The following relevant ClinVar are (likely) pathogenic: {skipped_exon_ClinVar.ids}"
        if not NMD_affected_exons:
            affected_exon = get_affected_exon(
                ref_transcript, transcript, variant, diff_len
            )[0]
        else:
            affected_exon = NMD_affected_exons[0]
        if is_NMD and path_disease_irrelevant_exons is not None:
            is_affected_exon_disease_relevant = check_exon_disease_relevant(
                path_disease_irrelevant_exons, NMD_affected_exons
            )
        else:
            is_affected_exon_disease_relevant = True
        is_reading_frame_preserved, _ = assess_reading_frame_preservation(diff_len)
        diff_len_protein_percent, ptc = calculate_prot_len_diff(
            ref_transcript, var_seq, diff_len
        )
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
            ptc=ptc,
            are_exons_skipped=are_exons_skipped,
            coding_exon_skipped=coding_exon_skipped,
            start_codon_exon_skipped=start_codon_exon_skipped,
            len_change_in_repetitive_region=len_change_in_repetitive_region,
            is_NMD=is_NMD,
            affected_exon=affected_exon,
            is_truncated_region_disease_relevant=is_truncated_region_disease_relevant,
            is_affected_exon_disease_relevant=is_affected_exon_disease_relevant,
            comment_truncated_region=comment_truncated_region_disease_relevant,
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
        try:
            ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
        except ValueError:
            raise Pyensembl_transcript_not_found
        try:
            ref_transcript.coding_sequence
        except ValueError or AttributeError:
            raise Pyensembl_no_coding_sequence
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
