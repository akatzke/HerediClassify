#!/usr/bin/env python3

import pathlib
import json

import pandas as pd

from cyvcf2 import VCF
from typing import Optional

from variant_classification.clinvar_utils import convert_vcf_gen_to_df


def convert_vcf_to_json(path_vcf: pathlib.Path) -> None:
    """
    For every entry in a VCF create a json file for input in classification tool
    """
    vcf = VCF(path_vcf)
    examples = convert_vcf_gen_to_df(vcf)
    for _, example in examples.iterrows():
        json_dict = create_json_dict_from_vcf(example)
        save_example_dict(json_dict, path_vcf)


def create_json_dict_from_vcf(data: pd.Series) -> dict:
    """
    Create and write json form vcf style pd.Series
    """
    json_dict = {}
    cons, gene, var_type = process_consequence(data.consequences)
    json_dict["chr"] = data.chrom
    json_dict["pos"] = int(data.pos)
    json_dict["gene"] = gene
    json_dict["ref"] = data.ref
    json_dict["alt"] = data.alt
    json_dict["variant_type"] = var_type
    json_dict["variant_effect"] = cons

    # Splicing prediction
    if type(data.spliceai_max_delta) is str:
        all_scores = []
        for spliceai_score in data.spliceai_max_delta.split(","):
            if spliceai_score != ".":
                all_scores.append(float(spliceai_score))
        if all_scores:
            json_dict["splicing_prediction_tools"] = {"SpliceAI": max(all_scores)}

    # Pathogenicity prediction
    pathogenicity_prediction = {}
    try:
        if isinstance(data.revel, str):
            revel_by_transcript = data.revel.split("|")
            revel_max_val = max([float(x.split("$")[1]) for x in revel_by_transcript])
            pathogenicity_prediction["REVEL"] = revel_max_val
        if isinstance(data.bayesdel, str):
            try:
                pathogenicity_prediction["BayesDel"] = float(data.bayesdel)
            except ValueError:
                scores_bayesdel = data.bayesdel.split("%26")
                float_scores = []
                for score in scores_bayesdel:
                    try:
                        float_scores.append(float(score))
                    except ValueError:
                        continue
                if float:
                    pathogenicity_prediction["BayesDel"] = max(float_scores)
        if len(pathogenicity_prediction) > 0:
            json_dict["pathogenicity_prediction_tools"] = pathogenicity_prediction
    except Exception:
        json_dict["pathogenicity_prediction_tools"] = pathogenicity_prediction

    # gnomAD
    gnomad_scores = {}
    try:
        if all([type(x) is str for x in [data.gnomad_af, data.gnomad_ac]]):
            gnomad_scores["AF"] = float(data.gnomad_af)
            gnomad_scores["AC"] = int(data.gnomad_ac)
            if all(
                [
                    type(x) is str
                    for x in [
                        data.gnomad_popmax,
                        data.gnomad_popmax_AC,
                        data.gnomad_popmax_AF,
                    ]
                ]
            ):
                gnomad_scores["subpopulation"] = data.gnomad_popmax
                gnomad_scores["popmax_AC"] = int(data.gnomad_popmax_AC)
                gnomad_scores["popmax_AF"] = float(data.gnomad_popmax_AF)
                try:
                    gnomad_scores["faf_popmax_AF"] = float(data.faf95_popmax)
                except AttributeError:
                    gnomad_scores["faf_popmax_AF"] = 0
            else:  # if the data is missing popmax values simply use the standard af and ac
                gnomad_scores["subpopulation"] = "ALL"
                gnomad_scores["popmax_AC"] = int(data.gnomad_ac)
                gnomad_scores["popmax_AF"] = float(data.gnomad_af)
                gnomad_scores["faf_popmax_AF"] = float(data.gnomad_af)
    except AttributeError:
        gnomad_scores = {}
    if len(gnomad_scores) > 0:
        json_dict["gnomAD"] = gnomad_scores

    # FLOSSIES
    flossies_score = {}
    try:
        if type(data.flossies_num_afr) is str:
            flossies_score["AFR"] = int(data.flossies_num_afr)
        if type(data.flossies_num_eur) is str:
            flossies_score["EUR"] = int(data.flossies_num_eur)
        if len(flossies_score) > 0:
            json_dict["FLOSSIES"] = flossies_score
    except AttributeError:
        flossies_score = {}

    # Cancer hotspots
    cancer_hotspots = {}
    try:
        if all(
            [type(x) is str for x in [data.cancerhotspots_af, data.cancerhotspots_ac]]
        ):
            cancer_hotspots["AF"] = float(data.cancerhotspots_af)
            cancer_hotspots["AC"] = int(data.cancerhotspots_ac)
        if len(cancer_hotspots) > 0:
            json_dict["cancer_hotspots"] = cancer_hotspots
    except AttributeError:
        json_dict["cancer_hotspots"] = cancer_hotspots

    # Cold spot
    json_dict["cold_spot"] = False

    # Functional and splice assay data
    try:
        if type(data.assays) is str:
            functional_assay, splice_assay = get_assay_data(data.assays)
            if functional_assay:
                json_dict["functional_data"] = functional_assay
            if splice_assay:
                json_dict["mRNA_analysis"] = splice_assay
    except AttributeError:
        return json_dict

    return json_dict


def process_consequence(cons: str) -> tuple[list[dict], str, list]:
    """
    Convert consequence to dictionary compatible with input
    Additionally returns gene affected by variant
    """
    cons_single = cons.split("&")
    cons_single_list = [con.split("|") for con in cons_single]
    keys = [
        "transcript",
        "hgvs_c",
        "hgvs_p",
        "variant_type",
        "impact",
        "exon",
        "intron",
        "gene",
        "protein_domain",
        "GENCODE_basics",
        "MANE_select",
        ".",
        "ensembl_canonical",
        "transcript_type",
        "length",
    ]
    cons_dict_list = []
    affected_genes = []
    for entry in cons_single_list:
        cons_dict = dict(zip(keys, entry))
        cons_dict_list.append(cons_dict)
        affected_genes.append(cons_dict["gene"])
    GOI = [
        "BRCA1",
        "BRCA2",
        "ATM",
        "CDH1",
        "PALB2",
        "PTEN",
        "TP53",
        "MSH6",
        "BARD1",
        "PMS2",
        "RAD51D",
        "PALB2",
        "MLH1",
        "APC",
        "RAD51B",
    ]
    try:
        gene = [gene for gene in GOI if gene in affected_genes][0]
    except IndexError:
        gene = max(affected_genes, key=affected_genes.count)
    gene_transcript_list = []
    selected_keys = [
        "transcript",
        "hgvs_c",
        "hgvs_p",
        "variant_type",
        "exon",
        "intron",
    ]
    for entry in cons_dict_list:
        if entry["gene"] == gene:
            select_transcript = {key: entry[key] for key in selected_keys}
            gene_transcript_list.append(select_transcript)
    reformatted_dict = reformat_consequence(gene_transcript_list)
    if not reformatted_dict:
        var_type = get_vartype_from_consequence(gene_transcript_list)
    else:
        var_type = get_vartype_from_consequence(reformatted_dict)
    return reformatted_dict, gene, var_type


def reformat_consequence(cons_list: list[dict]) -> list[dict]:
    """
    Reformat variabels in consequence dictionaries
    """
    out_cons_list = []
    for entry in cons_list:
        if entry["transcript"][0:4] == "ENST":
            if entry["hgvs_c"] == "None":
                continue
            entry["hgvs_c"] = entry["hgvs_c"].replace("%2B", "+")
            try:
                entry["exon"] = int(entry["exon"])
            except ValueError:
                del entry["exon"]
            try:
                entry["intron"] = int(entry["intron"])
            except ValueError:
                del entry["intron"]
            if entry["hgvs_p"] == "None":
                del entry["hgvs_p"]
            out_cons_list.append(entry)
            entry["variant_type"] = entry["variant_type"].split("_%26_")
    return out_cons_list


def get_vartype_from_consequence(cons_list_dict: list[dict]) -> list:
    """
    From the list of consequence, get all unique variant types
    """
    var_types = set()
    for cons_dict in cons_list_dict:
        if isinstance(cons_dict["variant_type"], str):
            var_type = [cons_dict["variant_type"]]
        else:
            var_type = cons_dict["variant_type"]
        var_types.update(var_type)
    return list(var_types)


def save_example_dict(json_dict: dict, in_path: pathlib.Path) -> None:
    """
    Save json dict to file_name created from gene and variant type
    """
    gene = json_dict["gene"]
    var_type = select_relevant_var_type(json_dict["variant_type"])
    file_name = gene + "_" + var_type + ".json"
    out_path = in_path.parent / file_name
    while out_path.exists():
        file_name = file_name.split(".json")[0] + "_2" + ".json"
        out_path = in_path.parent / file_name
    with open(out_path, "w") as f:
        json.dump(json_dict, f)


def select_relevant_var_type(var_types: list[str]) -> str:
    """
    From a list of variant types select relevant variant type
    """
    relevant_variant_types = [
        "stop_gained",
        "stop_lost",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "synonymous_variant",
        "start_lost",
        "splice_donor_variant",
        "splice_donor",
        "splice_acceptor_variant",
        "splice_acceptor",
    ]
    rel_var_type = [
        var_type for var_type in var_types if var_type in relevant_variant_types
    ]
    if rel_var_type:
        return rel_var_type[0]
    relevant_variant_types_bk = [
        "upstream_gene_variant",
        "downstream_gene_variant",
        "intron_variant",
        "3_prime_UTR_variant",
    ]
    rel_var_type_bk = [
        var_type for var_type in var_types if var_type in relevant_variant_types_bk
    ]
    if not rel_var_type_bk:
        print("This is the original var_type")
        print(var_types)
        return var_types[0]
    return rel_var_type_bk[0]


def get_assay_data(assay: str) -> Optional[tuple[list[dict], list[dict]]]:
    """
    Get assay data
    """
    if assay is None:
        return None
    if "&" in assay:
        assays = assay.split("&")
    else:
        assays = [assay]
    func_assay = []
    splice_assay = []
    for n in assays:
        assay_parts = n.split("|")
        if assay_parts[0] == "splicing":
            splice_dict = create_splice_dict(assay_parts[-1])
            splice_assay.append(splice_dict)
        elif assay_parts[0] == "functional":
            func_dict = create_func_dict(assay_parts[-1])
            func_assay.append(func_dict)
        else:
            raise ValueError
    return func_assay, splice_assay


def create_splice_dict(assay_info: str) -> dict:
    """
    From entry in form create results splicing assay
    """
    assay_info_dict = create_assay_info_dict(assay_info)
    out_dict = {}
    out_dict["minigene"] = assay_info_dict.get("Patient_RNA", False)
    out_dict["patient_rna"] = assay_info_dict.get("Minigene", False)
    out_dict["allelic"] = assay_info_dict.get("Allele-Specific", "False")
    out_dict["quantification"] = assay_info_dict.get(
        "Percent_aberrant_transcript", None
    )
    return out_dict


def create_func_dict(assay_info: str) -> dict:
    """
    From entry in form create results functional assay
    """
    assay_info_dict = create_assay_info_dict(assay_info)
    result_assay = assay_info_dict.get("Functional_category", None)
    if result_assay is None:
        print(
            f"There is an issue with the functional assay. The following results is given: {result_assay}."
        )
        out_dict = {"benign": False, "pathogenic": False}
    elif result_assay == "benign":
        out_dict = {"benign": True, "pathogenic": False}
    elif result_assay == "pathogenic":
        out_dict = {"benign": False, "pathogenic": True}
    elif result_assay == "ambigous":
        out_dict = {"benign": False, "pathogenic": False}
    else:
        print(
            f"There is an issue with the functional assay. The following results is given: {result_assay}."
        )
        out_dict = {"benign": False, "pathogenic": False}
    return out_dict


def create_assay_info_dict(assay_info: str) -> dict:
    assay_infos = assay_info.split("$")
    assay_info_dict = {
        entry.split("+")[0]: entry.split("+")[1] for entry in assay_infos
    }
    return assay_info_dict
