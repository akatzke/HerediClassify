#!/usr/bin/env python3

import pathlib

from cyvcf2 import VCF
import pyensembl
import pandas as pd

from clinvar_utils import (
    ClinVar,
    ClinVar_Type,
    convert_vcf_gen_to_df,
    filter_gene,
    create_ClinVar,
    summarise_ClinVars,
)
from variant import VariantInfo


def check_clinvar_start_alt_start(
    ref_transcript: pyensembl.transcript.Transcript,
    variant_info: VariantInfo,
    alt_start_codon: list[int],
    path_clinvar: pathlib.Path,
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
        variant_info, region_start, region_stop, path_clinvar
    )
    return ClinVar_start_alt_start


def check_clinvar_inframe_variant(
    variant: VariantInfo,
    path_clinvar_snv: pathlib.Path,
    path_clinvar_indel: pathlib.Path,
) -> ClinVar:
    start = variant.genomic_start
    stop = variant.genomic_end
    if start == stop:
        ClinVar_empty = create_ClinVar(pd.DataFrame(), ClinVar_Type.REGION)
        return ClinVar_empty
    clinvar_region_snv_df = get_clinvar_region_df(
        variant, start, stop, path_clinvar_snv
    )
    clinvar_region_indel_df = get_clinvar_region_df(
        variant, start, stop, path_clinvar_indel
    )
    clinvar_region_df = pd.concat([clinvar_region_snv_df, clinvar_region_indel_df])
    if not clinvar_region_df.empty:
        clinvar_region_df_lof = filter_lof_variants(clinvar_region_df)
        ClinVar_truncated_region = create_ClinVar(
            clinvar_region_df_lof, ClinVar_Type.REGION
        )
    else:
        ClinVar_truncated_region = create_ClinVar(pd.DataFrame(), ClinVar_Type.REGION)
    return ClinVar_truncated_region


def check_clinvar_truncated_region(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    path_clinvar_snv: pathlib.Path,
    path_clinvar_indel: pathlib.Path,
) -> ClinVar:
    """
    For exonic pvs1 variant coordinate appropriate clinvar entries
    """
    start, stop = define_range_truncation(ref_transcript, variant)
    clinvar_region_snv_df = get_clinvar_region_df(
        variant, start, stop, path_clinvar_snv
    )
    clinvar_region_indel_df = get_clinvar_region_df(
        variant, start, stop, path_clinvar_indel
    )
    clinvar_region_df = pd.concat([clinvar_region_snv_df, clinvar_region_indel_df])
    if not clinvar_region_df.empty:
        clinvar_region_df_lof = filter_lof_variants(clinvar_region_df)
        ClinVar_truncated_region = create_ClinVar(
            clinvar_region_df_lof, ClinVar_Type.REGION
        )
    else:
        ClinVar_truncated_region = create_ClinVar(pd.DataFrame(), ClinVar_Type.REGION)
    return ClinVar_truncated_region


def check_clinvar_NMD_exon(
    variant: VariantInfo,
    NMD_affected_exons: list[dict],
    path_clinvar_snv: pathlib.Path,
    path_clinvar_indel: pathlib.Path,
) -> ClinVar:
    """
    Check if exon contains any pathogenic variants
    """
    if NMD_affected_exons:
        ClinVar_exons = []
        for exon in NMD_affected_exons:
            clinvar_exon_snv_df = get_clinvar_region_df(
                variant, exon["exon_start"], exon["exon_end"], path_clinvar_snv
            )
            clinvar_exon_indel_df = get_clinvar_region_df(
                variant, exon["exon_start"], exon["exon_end"], path_clinvar_indel
            )
            clinvar_exon_df = pd.concat([clinvar_exon_snv_df, clinvar_exon_indel_df])
            if not clinvar_exon_df.empty:
                clinvar_exon_df_lof = filter_lof_variants(clinvar_exon_df)
                ClinVar_exon = create_ClinVar(clinvar_exon_df_lof, ClinVar_Type.REGION)
            else:
                ClinVar_exon = create_ClinVar(pd.DataFrame(), ClinVar_Type.REGION)
            ClinVar_exons.append(ClinVar_exon)
        ClinVar_exon_summary = summarise_ClinVars(
            ClinVar_exons, type=ClinVar_Type.REGION
        )
        return ClinVar_exon_summary
    else:
        ClinVar_exon = ClinVar(
            pathogenic=False, type=ClinVar_Type.REGION, highest_classification=None
        )
        return ClinVar_exon


def check_clinvar_region(
    variant_info: VariantInfo, start: int, end: int, path_clinvar: pathlib.Path
) -> ClinVar:
    """
    Get ClinVar entries in region
    """
    clinvar_region_filter = get_clinvar_region_df(
        variant_info, start, end, path_clinvar
    )
    ClinVar_region = create_ClinVar(clinvar_region_filter, ClinVar_Type.REGION)
    return ClinVar_region


def get_clinvar_region_df(
    variant_info: VariantInfo, start: int, end: int, path_clinvar: pathlib.Path
) -> pd.DataFrame:
    """
    Get ClinVar dataframe for region
    Filtered for gene of interest
    """
    clinvar = VCF(path_clinvar)
    clinvar_region = clinvar(f"{variant_info.chr}:{start}-{end}")
    clinvar_region_df = convert_vcf_gen_to_df(clinvar_region)
    if clinvar_region_df.empty:
        return pd.DataFrame()
    clinvar_region_filter = filter_gene(clinvar_region_df, variant_info.gene_name)
    return clinvar_region_filter


def filter_lof_variants(clinvar_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for loss of function variants in MC column
    """
    drop_na = clinvar_df.dropna(subset=["MC"])
    lof_variants = drop_na[
        drop_na.MC.str.contains("frameshift") | drop_na.MC.str.contains("nonsense")
    ]
    return lof_variants


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
