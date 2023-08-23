#!/usr/bin/env python3

from variant import TranscriptInfo, Variant
from typing import List
from recycled_code_genotoscope import assess_exon_skipping


def is_transcript_disease_relevant(
    variant: Variant, disease_relevant_transcripts: List[str]
) -> None:
    for transcript in variant.transcript_info:
        if transcript.transcript_id in disease_relevant_transcripts:
            transcript.transcript_disease_relevant = True
        else:
            transcript.transcript_disease_relevant = False


def is_len_change_in_repetitive_region(
    transcript: TranscriptInfo, repetitive_regions
) -> None:
    if (
        transcript.var_start in repetitive_regions
        and transcript.var_stop in repetitive_regions
    ):
        transcript.len_change_in_repetitive_region = True
    else:
        transcript.len_change_in_repetitive_region = False


def assess_effect_splice_variants(variant: Variant) -> None:
    for transcript in variant.transcript_info:
        if any(var_type in transcript.var_type for var_type in ["splicing_phenotypes"]):
            assess_exon_skipping(transcript, variant.variant_info)
        else:
            continue
