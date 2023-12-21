#!/usr/bin/env python3

import logging

from variant import Variant

logger = logging.getLogger("GenOtoScope_Classify.check_disease_relenvat_transcript")


def check_disease_relevant_transcript(variant: Variant, config: dict) -> Variant:
    """
    Check if a disease relevant transcript was defined in the configuration
    If so remove all other transcripts
    If not return same object
    """
    try:
        disease_relevant_transcripts = config["disease_relevant_transcripts"]
    except KeyError:
        logger.debug(
            f"No disease relevant transcript defined in configuration {config['name']} for {variant.variant_info.gene_name}."
        )
        return variant
    disease_relevant_transcript_info = []
    for transcript in variant.transcript_info:
        if transcript.transcript_id in [
            transcript["name"] for transcript in disease_relevant_transcripts
        ]:
            disease_relevant_transcript_info.append(transcript)
    if not disease_relevant_transcript_info:
        logger.debug(
            f"The disease relevant transcript {disease_relevant_transcripts} was not found in list of transcripts for the variant {[transcript.transcript_id for transcript in variant.transcript_info]}. Use complete list of transcripts instead."
        )
        return variant
    variant.transcript_info = disease_relevant_transcript_info
    return variant
