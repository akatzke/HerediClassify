#!/usr/bin/env python3

import pathlib
from typing import Optional

import pyensembl

import pandas as pd

from pybedtools import BedTool, Interval
from variant import TranscriptInfo, VariantInfo


def check_bed_intersect_start_loss(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    alt_start_codon: list[int],
    path_rep_uniprot: pathlib.Path,
) -> tuple[bool, str]:
    """
    Assess if prot_len_change caused by alternative start codon is in repetitive region
    """
    ref_start_codon = ref_transcript.start_codon_positions
    if ref_transcript.strand == "-":
        gen_start = min(alt_start_codon)
        gen_end = max(ref_start_codon)
    else:
        gen_start = min(ref_start_codon)
        gen_end = max(alt_start_codon)
    prot_len_in_repetitive_region, comment = check_intersection_with_bed(
        variant, gen_start, gen_end, ref_transcript, path_rep_uniprot
    )
    return prot_len_in_repetitive_region, comment


def check_intersection_with_bed(
    variant: VariantInfo,
    gen_start: int,
    gen_end: int,
    ref_transcript: pyensembl.transcript.Transcript,
    path_bed: pathlib.Path,
) -> tuple[bool, str]:
    """
    Check if variant overlaps UniProt annotated repetitive region
    """
    variant_interval = BedTool(
        create_bed_line(variant, gen_start, gen_end, ref_transcript.strand),
        from_string=True,
    )[0]
    bed = BedTool(path_bed).sort()
    annotation_hits = bed.all_hits(variant_interval, same_strand=True)
    if len(annotation_hits) > 0:
        comment = create_comment_from_bed_file(annotation_hits, path_bed)
        return True, comment
    return False, ""


def create_bed_line(
    variant: VariantInfo, gen_start: int, gen_end: int, transcript_strand: str
) -> str:
    """
    Create bed line to represent variant info
    Credits: http://daler.github.io/pybedtools/intervals.html
    """
    bed_line = " ".join(
        [
            "chr" + str(variant.chr),
            str(gen_start),
            str(gen_end),
            variant.gene_name,
            ".",
            transcript_strand,
        ]
    )
    return bed_line


def create_comment_from_bed_file(hits: list[Interval], path_bed: pathlib.Path) -> str:
    """
    From intersections with bed file, create a comment
    """
    bed = pd.read_csv(path_bed, sep="\t")
    hits_df = pd.DataFrame()
    for hit in hits:
        chr = hit[0]
        start = int(hit[1])
        end = int(hit[2])
        bed_entries = bed[
            (bed["#chr"] == chr) & (bed["start"] == start) & (bed["end"] == end)
        ]
        hits_df = pd.concat([hits_df, bed_entries])
    if "domain_name" in bed.columns:
        hits_df["comment"] = hits_df["domain_name"] + " in " + hits_df["gene"]
        comment = " ".join(hits_df.comment)
    elif "comment" in bed.columns:
        comment = " ".join(hits_df.comment)
    else:
        comment = ""
    return comment


def select_mane_transcript(
    transcripts: list[TranscriptInfo], mane_path: pathlib.Path
) -> Optional[TranscriptInfo]:
    """
    From a list of transcripts select the MANE transcript
    """
    mane_transcripts_df = pd.read_csv(mane_path, sep="\t")
    mane_transcripts = mane_transcripts_df.transcript.dropna()
    for transcript in transcripts:
        if any(mane_transcripts.str.contains(transcript.transcript_id)):
            return transcript
    return None
