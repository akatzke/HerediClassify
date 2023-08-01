#!/usr/bin/env python3

import logging

logger = logging.getLogger("GenOtoScope_Classify.PVS1.exon_skipping")


def examine_alternative_start_codon(ref_transcript) -> bool:
    """
    Check on the other transcripts of the same gene, if they use an alternative start codon
    """

    logger.debug(
        "Check if alternate transcripts, of the same gene containing the variant, use alternative start codon positions"
    )
    exists_alternative_start_codon = True
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

    # if all alternate transcripts have the same start codon position,
    # use the same chromosomal position for start codon
    # print(alternate_transcripts_unique_start_chr_pos)
    if len(
        alternate_transcripts_unique_start_chr_pos
    ) == alternate_transcripts_unique_start_chr_pos.count(var_start_codon_chr_pos):
        exists_alternative_start_codon = False
    else:
        exists_alternative_start_codon = True
    return exists_alternative_start_codon
