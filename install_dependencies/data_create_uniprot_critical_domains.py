#!/usr/bin/env python3

import pathlib
import argparse
import cyvcf2
import pandas as pd
import numpy as np

from cyvcf2 import VCF

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


def create_uniprot_critical_domains(
    uniprot_path: pathlib.Path, clinvar_path: pathlib.Path, out_path: pathlib.Path
):
    """
    From UniProt and ClinVar create disease relevant transcript
    """
    uniprot = pd.read_csv(uniprot_path, sep="\t")
    clinvar = VCF(clinvar_path)
    uniprot.rename(columns={"#chrom": "chrom"}, inplace=True)
    uniprot_out = uniprot.assign(
        count_lp_p=np.nan,
        count_lb_b=np.nan,
        count_vus=np.nan,
        pathogenic_probability=np.nan,
    )
    uniprot_patho_score = uniprot_out.apply(
        calculate_pathogenicity_score, clinvar=clinvar, axis=1
    )


def calculate_pathogenicity_score(
    domain: pd.Series, clinvar: cyvcf2.cyvcf2.VCF
) -> pd.Series:
    """
    Calculate pathogenicity score for a domain
    """
    overlapping_clinvar = get_clinvar_entries(clinvar, domain)
    count_by_classification = get_count_by_classification(overlapping_clinvar)
    out = create_output_entry(domain, count_by_classification)
    return out


def get_clinvar_entries(clinvar: cyvcf2.cyvcf2.VCF, domain: pd.Series) -> pd.DataFrame:
    """
    From ClinVar get entries located in domain
    """
    clinvar_region = clinvar(
        f"{domain.chrom.split('chr')[1]}:{domain.chromStart}-{domain.chromEnd}"
    )
    clinvar_region_df = convert_vcf_gen_to_df(clinvar_region)
    return clinvar_region_df


def get_count_by_classification(clinvars: pd.DataFrame) -> dict[str:int]:
    quality_two_star = [
        "reviewed_by_expert_panel",
        "criteria_provided,_multiple_submitters,_no_conflicts",
    ]
    clinvars_quality = clinvars[clinvars.CLNREVSTAT.isin(quality_two_star)]


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
