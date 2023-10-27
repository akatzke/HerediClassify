#!/usr/bin/env python3

import logging
import pathlib
from collections.abc import Iterable

import hgvs.parser
import pyensembl
from cyvcf2 import VCF

from refactoring.variant import VariantInfo, TranscriptInfo
from refactoring.var_type import VARTYPE_GROUPS
from refactoring.clinvar_utils import (
    ClinVar,
    CLINVAR_TYPE,
    convert_vcf_gen_to_df,
    create_ClinVar,
    get_affected_transcript,
)
from refactoring.genotoscope_exon_skipping import (
    parse_variant_intron_pos,
    find_exon_by_ref_pos,
)

hgvs_parser = hgvs.parser.Parser()

logger = logging.getLogger("GenOtoScope_Classify.clinvar.splicing")

# path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730.vcf.gz")


def check_clinvar_splicing(
    variant: VariantInfo,
    transcripts: Iterable[TranscriptInfo],
    path_clinvar: pathlib.Path,
) -> tuple[ClinVar, ClinVar]:
    """
    Check ClinVar for entries supporting pathogenicity of splice site
    """
    clinvar = VCF(path_clinvar)
    ### Check ClinVar for pathogenic variants with same nucleotide change
    clinvar_same_pos = clinvar(
        f"{variant.chr}:{variant.genomic_start}-{variant.genomic_end}"
    )
    clinvar_same_pos_df = convert_vcf_gen_to_df(clinvar_same_pos)
    ClinVar_same_pos = create_ClinVar(clinvar_same_pos_df, CLINVAR_TYPE.SAME_NUCLEOTIDE)
    ### Check ClinVar for pathogenic variant in same / closest splice site
    affected_transcript = get_affected_transcript(transcripts, VARTYPE_GROUPS.INTRONIC)
    (start_splice_site, end_splice_site) = find_corresponding_splice_site(
        affected_transcript, variant
    )
    clinvar_splice_site = clinvar(
        f"{variant.chr}:{start_splice_site}-{end_splice_site}"
    )
    clinvar_splice_site_df = convert_vcf_gen_to_df(clinvar_splice_site)
    ClinVar_splice_site = create_ClinVar(
        clinvar_splice_site_df, CLINVAR_TYPE.SAME_SPLICE_SITE
    )
    return (ClinVar_same_pos, ClinVar_splice_site)


def find_corresponding_splice_site(
    transcript: TranscriptInfo, variant: VariantInfo
) -> tuple[int, int]:
    """
    Reconstruct splice site
    Splice site is defined as +/- 1,2 as only for these locations varinat is clearly defines as a splice variant
    """
    ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
        transcript.transcript_id
    )
    if "+" in str(transcript.var_hgvs) or "-" in str(transcript.var_hgvs):
        splice_site_start, splice_site_stop = get_splice_site_for_intronic_variant(
            variant, transcript, ref_transcript
        )
    else:
        splice_site_start, splice_site_stop = get_splice_site_for_exonic_variant(
            variant, ref_transcript
        )
    return (splice_site_start, splice_site_stop)


def get_splice_site_for_intronic_variant(
    variant: VariantInfo,
    transcript: TranscriptInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple[int, int]:
    """
    Get splice site for intronic variant
    """
    (
        split_symbol,
        distance_to_splice_site,
        direction_to_splice_site,
    ) = parse_variant_intron_pos(transcript.var_hgvs)
    distance_to_splice_site_start = distance_to_splice_site - 1
    distance_to_splice_site_stop = distance_to_splice_site - 2
    if ref_transcript.strand == "-":
        # Reverse direction_to_splice_site, to point into genomic direction of splice site
        direction_to_splice_site = direction_to_splice_site * -1
    splice_site_start = (
        variant.genomic_start + direction_to_splice_site * distance_to_splice_site_start
    )
    splice_site_stop = (
        variant.genomic_start + direction_to_splice_site * distance_to_splice_site_stop
    )
    if splice_site_start > splice_site_stop:
        return splice_site_stop, splice_site_start
    else:
        return splice_site_start, splice_site_stop


def get_splice_site_for_exonic_variant(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple[int, int]:
    """
    Get splice site for exonic variant
    """
    exon_index, pos_in_exon = find_exon_by_ref_pos(
        ref_transcript, variant.genomic_start, is_genomic=True
    )
    affected_exon = ref_transcript.exons[exon_index]
    if abs(variant.genomic_start - affected_exon.start) < abs(
        variant.genomic_start - affected_exon.end
    ):
        splice_site_start = affected_exon.start - 2
        splice_site_end = affected_exon.start - 1
    else:
        splice_site_start = affected_exon.end + 1
        splice_site_end = affected_exon.end + 2
    return splice_site_start, splice_site_end
