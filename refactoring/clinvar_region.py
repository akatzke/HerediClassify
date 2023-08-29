#!/usr/bin/env python3

import pathlib

from cyvcf2 import VCF
import pyensembl

from refactoring.clinvar_utils import convert_vcf_gen_to_df, filter_gene, create_ClinVar
from refactoring.variant import VariantInfo
from refactoring.variant_annotate import ClinVar


def check_clinvar_start_alt_start(
    ref_transcript: pyensembl.transcript.Transcript,
    variant_info: VariantInfo,
    alt_start_codon: list[int],
) -> ClinVar:
    """
    Get ClinVar entries between start codon and closest alternative start_codon
    """
    path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730.vcf.gz")
    ref_start_codon = ref_transcript.start_codon_positions
    if ref_transcript.strand == "-":
        region_start = alt_start_codon[2] + 1
        region_stop = ref_start_codon[2]
    else:
        region_start = ref_start_codon[0]
        region_stop = alt_start_codon[0] - 1
    ClinVar_start_alt_start = check_clinvar_region(
        variant_info, region_start, region_stop, path_clinvar
    )
    return ClinVar_start_alt_start


def check_clinvar_NMD_exon(variant: VariantInfo, NMD_affected_exon: dict) -> ClinVar:
    """
    Check if exon contains any pathogenic variants
    """
    path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730.vcf.gz")
    ClinVar_exon = check_clinvar_region(
        variant, NMD_affected_exon["start"], NMD_affected_exon["end"], path_clinvar
    )
    return ClinVar_exon


def check_clinvar_region(
    variant_info: VariantInfo, start: int, end: int, path_clinvar: str
) -> ClinVar:
    """
    Get ClinVar entries in region
    """
    clinvar = VCF(path_clinvar)
    clinvar_region = clinvar(f"{variant_info.chr}:{start}-{end}")
    clinvar_region_df = convert_vcf_gen_to_df(clinvar_region)
    clinvar_region_filter = filter_gene(clinvar_region_df, variant_info.gene_name)
    ClinVar_region = create_ClinVar(clinvar_region_filter, "region")
    return ClinVar_region
