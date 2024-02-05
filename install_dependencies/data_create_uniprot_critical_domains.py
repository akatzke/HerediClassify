#!/usr/bin/env python3

import pathlib
import argparse
import cyvcf2
import pandas as pd
import numpy as np

from cyvcf2 import VCF
import pyensembl

from variant_classification.clinvar_utils import convert_vcf_gen_to_df


parser = argparse = argparse.ArgumentParser()

parser.add_argument(
    "-u",
    "--uniprot",
    default="",
    help="path to uniprot domain file",
    type=str,
)
parser.add_argument(
    "-c", "--clinvar", default="", help="path to unfiltered clinvar file", type=str
)
parser.add_argument(
    "-o",
    "--output",
    default="",
    help="output file path. If not given will default to standard output.",
)
args = parser.parse_args()

ensembl = pyensembl.EnsemblRelease(110)


def create_uniprot_critical_domains(
    uniprot_path: pathlib.Path, clinvar_path: pathlib.Path, out_path: pathlib.Path
):
    """
    From UniProt and ClinVar create disease relevant transcript
    """
    clinvar = VCF(clinvar_path)
    uniprot = pd.read_csv(uniprot_path, sep="\t")
    uniprot.rename(columns={"#chrom": "chrom"}, inplace=True)
    uniprot_chrom = uniprot[uniprot.chrom.str.match("chr[0-9]+$")]
    uniprot_out = uniprot_chrom.assign(
        count_lp_p=np.nan,
        count_lb_b=np.nan,
        count_vus=np.nan,
        pathogenicity_score=np.nan,
    )
    uniprot_patho_score = uniprot_out.apply(
        annotate_pathogenicity_score, clinvar=clinvar, axis=1
    )
    critical_domains = uniprot_patho_score[
        (uniprot_patho_score.pathogenicity_score >= 0.51)
        & (uniprot_patho_score.count_lp_p >= 5)
        & (uniprot_patho_score.count_lb_b == 0)
    ]
    return critical_domains


def annotate_pathogenicity_score(
    domain: pd.Series, clinvar: cyvcf2.cyvcf2.VCF
) -> pd.Series:
    """
    Calculate pathogenicity score for a domain
    """
    overlapping_clinvar = get_clinvar_entries(clinvar, domain)
    if overlapping_clinvar.empty:
        domain.pathogenicity_score = 0
        domain.count_lp_p = 0
        domain.count_lb_b = 0
        domain.count_vus = 0
        return domain
    count_by_class = get_count_by_classification(overlapping_clinvar)
    pathogenicity_score = calculate_pathogenicity_score(count_by_class)
    domain.pathogenicity_score = pathogenicity_score
    domain.count_lp_p = count_by_class["num_p_lp"]
    domain.count_lb_b = count_by_class["num_b_lb"]
    domain.count_vus = count_by_class["num_VUS"]
    return domain


def get_clinvar_entries(clinvar: cyvcf2.cyvcf2.VCF, domain: pd.Series) -> pd.DataFrame:
    """
    From ClinVar get entries located in domain
    """
    clinvar_region = clinvar(
        f"{domain.chrom.split('chr')[1]}:{domain.chromStart}-{domain.chromEnd}"
    )
    clinvar_region_df = convert_vcf_gen_to_df(clinvar_region)
    if clinvar_region_df.empty:
        return clinvar_region_df
    quality_two_star = [
        "reviewed_by_expert_panel",
        "criteria_provided,_multiple_submitters,_no_conflicts",
    ]
    clinvar_quality = clinvar_region_df[
        clinvar_region_df.CLNREVSTAT.isin(quality_two_star)
    ]
    if clinvar_quality.empty:
        return clinvar_quality
    clinvar_filtered = filter_by_strand(clinvar_quality, domain.strand)
    return clinvar_filtered


def filter_by_strand(clinvar_entries: pd.DataFrame, strand: str) -> pd.DataFrame:
    """
    Check for all clinvar entries if gene is on the same strand as the protein domain
    """
    # Format GENEINFO column to only contain gene name
    clinvar_entries["gene"] = clinvar_entries.GENEINFO.str.split(":", expand=True)[0]
    unique_genes = list(clinvar_entries.gene.unique())
    genes_to_remove = []
    for gene in unique_genes:
        try:
            ensembl_genes = ensembl.genes_by_name(gene)
        except ValueError:
            continue
        for ensembl_gene in ensembl_genes:
            if ensembl_gene.strand != strand:
                genes_to_remove.append(gene)
    filtered_clinvar = clinvar_entries[~clinvar_entries.gene.isin(genes_to_remove)]
    return filtered_clinvar


def get_count_by_classification(clinvar: pd.DataFrame) -> dict[str, int]:
    """
    Filter ClinVar for entreis
    """
    num_VUS = clinvar[clinvar.CLNSIG == "Uncertain_significance"].shape[0]
    num_b_lb = clinvar[clinvar.CLNSIG.str.contains("benign", case=False)].shape[0]
    num_p_lp = clinvar[clinvar.CLNSIG.str.contains("pathogenic", case=False)].shape[0]
    return {"num_VUS": num_VUS, "num_b_lb": num_b_lb, "num_p_lp": num_p_lp}


def calculate_pathogenicity_score(count_by_class: dict[str, int]) -> float:
    """
    Calcualte pathogenicity score
    """
    delta = 10 ** (-6)
    total = sum(
        [
            count_by_class["num_b_lb"] + delta,
            count_by_class["num_p_lp"] + delta,
            count_by_class["num_VUS"] + delta,
        ]
    )
    return (count_by_class["num_p_lp"] + delta) / total


def main():
    if args.uniprot == "":
        raise ValueError("No uniprot file provided.")
    uniprot_path = pathlib.Path(args.uniprot)
    if args.clinvar == "":
        raise ValueError("No clinvar file provided.")
    clinvar_path = pathlib.Path(args.clinvar)
    if args.output == "":
        file_name = f"{uniprot_path.stem}.bed"
        out_path = uniprot_path.parent / file_name
    else:
        out_path = pathlib.Path(args.output)
    create_uniprot_critical_domains(uniprot_path, clinvar_path, out_path)


if __name__ == "__main__":
    main()
