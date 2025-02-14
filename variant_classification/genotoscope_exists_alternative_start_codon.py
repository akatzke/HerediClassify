#!/usr/bin/env python3

import logging
import pyensembl

from variant import VariantInfo
from genotoscope_exon_skipping import (
    extract_codons,
    is_transcript_in_positive_strand,
    find_exon_by_ref_pos,
)

logger = logging.getLogger("GenOtoScope_Classify.PVS1.find_alternative_start_codon")


def assess_alternative_start_codon(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    var_coding_seq: str,
) -> tuple[bool, list[int], list[int]]:
    (
        alternative_start_codons_cDNA,
        alternative_start_codons_genomic,
    ) = find_alternative_start_codons(variant, ref_transcript, var_coding_seq)
    if alternative_start_codons_genomic == [0, 0, 0]:
        return False, [0, 0, 0], [0, 0, 0]
    elif alternative_start_codons_cDNA == [0, 1, 2]:
        return False, [0, 0, 0], [0, 0, 0]
    else:
        alternative_start_codons_genomic.sort()
        return True, alternative_start_codons_genomic, alternative_start_codons_cDNA


def find_alternative_start_codons(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    var_coding_seq: str,
) -> tuple[list[int], list[int]]:
    """
    First seach for alternative start codon in other reference transcripts
    If no start codon is found, search in sequence for other start codons
    First seach in the first 200 codons, then search in whole sequence
    Return genomic positions of start codon
    """
    downstream_codons = get_codons_downstream_start(var_coding_seq, 201)
    codon_position_start_codon = search_closest_start_codon(downstream_codons, variant)
    if 0 > codon_position_start_codon:
        all_downstream_codons = get_codons_downstream_start(
            var_coding_seq, len(var_coding_seq)
        )
        codon_position_start_codon = search_closest_start_codon(
            all_downstream_codons, variant
        )
    if codon_position_start_codon >= 0:
        alternative_start_codon = covnert_codon_position_to_genomic_position(
            ref_transcript, codon_position_start_codon
        )
        cdna_alternative_start_codon = get_nucleotide_positions_from_codon_position(
            codon_position_start_codon
        )
        return cdna_alternative_start_codon, alternative_start_codon
    else:
        return [0, 0, 0], [0, 0, 0]


def examine_start_codon_other_transcripts(ref_transcript) -> list[int]:
    """
    Check on the other transcripts of the same gene, if they use an alternative start codon
    Return the closest upstream start codon
    """

    logger.debug(
        "Check if alternate transcripts, of the same gene containing the variant, use alternative start codon positions"
    )
    ref_transcript_id = ref_transcript.id
    var_start_codon_chr_pos = ref_transcript.start_codon_positions
    logger.debug(
        f"variant transcript start codon chr positions: {var_start_codon_chr_pos}"
    )
    alternate_transcripts_unique_start_chr_pos = []
    for transcript in ref_transcript.gene.transcripts:
        if (
            transcript.id != ref_transcript_id
        ):  # if transcript is different that the one harbouring the variant
            if (
                transcript.contains_start_codon
            ):  # if start codon is annotated in the transcript
                if (
                    transcript.start_codon_positions
                    not in alternate_transcripts_unique_start_chr_pos
                ):
                    alternate_transcripts_unique_start_chr_pos.append(
                        transcript.start_codon_positions
                    )
    alternate_transcripts_unique_start_chr_pos.remove(var_start_codon_chr_pos)
    # Find closest start codon in alternative transcript
    if not alternate_transcripts_unique_start_chr_pos:
        return []
    elif is_transcript_in_positive_strand(ref_transcript):
        sorted_start_codons = sorted(alternate_transcripts_unique_start_chr_pos)
        for start_codon in sorted_start_codons:
            if min(start_codon) > min(ref_transcript.start_codon_positions):
                return start_codon
        return []
    else:
        sorted_start_codons = sorted(
            alternate_transcripts_unique_start_chr_pos, reverse=True
        )
        for start_codon in sorted_start_codons:
            if max(start_codon) < max(ref_transcript.start_codon_positions):
                return start_codon
        return []


def get_codons_downstream_start(
    var_coding_seq: str, downstream_length: int
) -> list[str]:
    """
    Get codons in downstream region from start codon of transcript
    """

    return extract_codons(var_coding_seq[0:downstream_length])


def search_closest_start_codon(codons: list[str], variant: VariantInfo) -> int:
    """
    Search for start codon in input codons
    if exists, return the closest to the start

    positive value = index of start codon in input list of codons,
    negative value = start codon is not contained in the input list of codons
    """

    try:
        return codons.index("ATG")
    except ValueError:
        logger.debug(
            f"Codons do not contain start codon\n=> variant position: {variant.to_string()}"
        )
        return -1


def covnert_codon_position_to_genomic_position(
    ref_transcript: pyensembl.transcript.Transcript, alt_start_codon_index: int
) -> list[int]:
    """
    Covert start codon index into list of genomic positions of start codon
    """

    # print("Extract pathogenic variants in upstream of the closest codon")
    is_genomic = False
    closest_start_codon_pos = (alt_start_codon_index) * 3
    (
        closest_start_codon_exon_index,
        closest_start_codon_exon_offset,
    ) = find_exon_by_ref_pos(ref_transcript, closest_start_codon_pos, is_genomic)
    # convert exon position in genomic position for variants retrieval
    closest_start_codon_pos = convert_exon_pos2genomic_pos(
        ref_transcript, closest_start_codon_exon_index, closest_start_codon_exon_offset
    )
    if is_transcript_in_positive_strand(ref_transcript):
        complete_start_codon = [
            closest_start_codon_pos - 2,
            closest_start_codon_pos - 1,
            closest_start_codon_pos,
        ]
    else:
        complete_start_codon = [
            closest_start_codon_pos,
            closest_start_codon_pos + 1,
            closest_start_codon_pos + 2,
        ]

    return complete_start_codon


def convert_exon_pos2genomic_pos(
    ref_transcript: pyensembl.transcript.Transcript,
    overlap_exon_index: int,
    overlap_exon_offset: int,
) -> int:
    """
    Convert overlap position in exonic coordinates to genomic position
    """

    overlap_exon = ref_transcript.exons[overlap_exon_index]

    if is_transcript_in_positive_strand(ref_transcript):
        # for positive strand, add to the exon start
        exon_start = overlap_exon.to_dict()["start"]
        overlap_genomic_pos = exon_start + overlap_exon_offset
    else:
        # for negative strand, subtract from the exon end
        exon_end = overlap_exon.to_dict()["end"]
        overlap_genomic_pos = exon_end - overlap_exon_offset
    return overlap_genomic_pos


def get_nucleotide_positions_from_codon_position(codon_pos: int) -> list[int]:
    """
    From the codon position get the position of all nucleotides in the codon
    """
    start_nuc = codon_pos * 3
    mid_nuc = start_nuc + 1
    end_nuc = start_nuc + 2
    return [start_nuc, mid_nuc, end_nuc]
