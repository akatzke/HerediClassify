#!/usr/bin/env python3

import pathlib

from cyvcf2 import VCF
import pyensembl

from refactoring.clinvar_utils import (
    ClinVar,
    CLINVAR_TYPE,
    convert_vcf_gen_to_df,
    filter_gene,
    create_ClinVar,
    summarise_ClinVars,
)
from refactoring.variant import VariantInfo


def check_clinvar_start_alt_start(
    ref_transcript: pyensembl.transcript.Transcript,
    variant_info: VariantInfo,
    alt_start_codon: list[int],
) -> ClinVar:
    """
    Get ClinVar entries between start codon and closest alternative start_codon
    """
    ref_start_codon = ref_transcript.start_codon_positions
    if ref_transcript.strand == "-":
        region_start = alt_start_codon[2] + 1
        region_stop = ref_start_codon[2]
    else:
        region_start = ref_start_codon[0]
        region_stop = alt_start_codon[0] - 1
    ClinVar_start_alt_start = check_clinvar_region(
        variant_info, region_start, region_stop
    )
    return ClinVar_start_alt_start


def check_clinvar_truncated_region(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> ClinVar:
    """
    For exonic pvs1 variant coordinate appropriate clinvar entries
    """
    start, stop = define_range_truncation(ref_transcript, variant)
    clinvar_exonic = check_clinvar_region(variant, start, stop)
    return clinvar_exonic


def check_clinvar_NMD_exon(
    variant: VariantInfo, NMD_affected_exons: list[dict]
) -> ClinVar:
    """
    Check if exon contains any pathogenic variants
    """
    if NMD_affected_exons:
        ClinVar_exons = []
        for exon in NMD_affected_exons:
            ClinVar_exon = check_clinvar_region(
                variant, exon["exon_start"], exon["exon_end"]
            )
            ClinVar_exons.append(ClinVar_exon)
        ClinVar_exon_summary = summarise_ClinVars(ClinVar_exons, type= CLINVAR_TYPE.REGION)
        return ClinVar_exon_summary
    else:
        ClinVar_exon = ClinVar(
            pathogenic=False, type= CLINVAR_TYPE.REGION, highest_classification=None
        )
        return ClinVar_exon


def check_clinvar_region(variant_info: VariantInfo, start: int, end: int) -> ClinVar:
    """
    Get ClinVar entries in region
    """
    path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730_snv.vcf.gz")
    clinvar = VCF(path_clinvar)
    clinvar_region = clinvar(f"{variant_info.chr}:{start}-{end}")
    clinvar_region_df = convert_vcf_gen_to_df(clinvar_region)
    clinvar_region_filter = filter_gene(clinvar_region_df, variant_info.gene_name)
    ClinVar_region = create_ClinVar(clinvar_region_filter, CLINVAR_TYPE.REGION)
    return ClinVar_region


def define_range_truncation(
    ref_transcript: pyensembl.transcript.Transcript, variant: VariantInfo
) -> tuple[int, int]:
    """
    Based on strand construct complete region affected by truncating variant
    """
    if ref_transcript.strand == "+":
        start = variant.genomic_start
        end = ref_transcript.stop_codon_positions[-1]
        return (start, end)
    else:
        end = variant.genomic_end
        start = ref_transcript.stop_codon_positions[0]
        return (start, end)


def define_range_truncation_exon(
    ref_transcript: pyensembl.transcript.Transcript, variant: VariantInfo
) -> tuple[int, int]:
    """
    Based on strand construct region in exon affected by truncation
    """
    if ref_transcript.strand == "+":
        start = variant.genomic_start
        end = ref_transcript.exon.stop
        return (start, end)
    else:
        end = variant.genomic_end
        start = ref_transcript.exon.stop
        return (start, end)
