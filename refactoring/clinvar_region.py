#!/usr/bin/env python3

import pathlib

from cyvcf2 import VCF
import pyensembl

from refactoring.clinvar_utils import (
    convert_vcf_gen_to_df,
    filter_gene,
    create_ClinVar,
    summarise_ClinVars,
)
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


def check_clinvar_NMD_exon(variant: VariantInfo, NMD_affected_exon: dict) -> ClinVar:
    """
    Check if exon contains any pathogenic variants
    """
    if NMD_affected_exon:
        ClinVar_exons = []
        for exon in NMD_affected_exon:
            ClinVar_exon = check_clinvar_region(
                variant, exon["exon_start"], exon["exon_end"]
            )
            ClinVar_exons.append(ClinVar_exon)
        ClinVar_exon_summary = summarise_ClinVars(ClinVar_exons, type="region")
        return ClinVar_exon_summary
    else:
        ClinVar_exon = ClinVar(
            pathogenic=False, type="region", highest_classification=None, ids=[]
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
    ClinVar_region = create_ClinVar(clinvar_region_filter, "region")
    return ClinVar_region
