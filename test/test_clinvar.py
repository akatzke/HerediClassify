#!/usr/bin/env python3

import pathlib

import test.paths as paths
from test.test_import_variant import create_json_string_from_variant

from variant_classification.clinvar_annot import annotate_clinvar
from variant_classification.clinvar_utils import ClinVar_Type
from variant_classification.load_config import load_config
from variant_classification.load_variant import load_variant


def test_missense_2():
    """
    Test ClinVar annotation for missense variant
    """
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "BRCA1_missense_variant.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    # Create needed paths
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_files = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_files / file_name
    clinvar_dict = annotate_clinvar(
        variant.variant_info, variant.transcript_info, path_clinvar
    )
    for clinvar_type, clinvar_object in clinvar_dict.items():
        if clinvar_type.value == ClinVar_Type.SAME_AA_CHANGE.value:
            clinvar_same_aa_change = clinvar_object
        elif clinvar_type.value == ClinVar_Type.DIFF_AA_CHANGE.value:
            clinvar_diff_aa_change = clinvar_object
        elif clinvar_type.value == ClinVar_Type.SAME_NUCLEOTIDE.value:
            clinvar_same_nucleotide = clinvar_object
        elif clinvar_type.value == ClinVar_Type.SAME_SPLICE_SITE.value:
            clinvar_same_splice_site = clinvar_object
        else:
            continue
    assert (
        clinvar_same_aa_change.pathogenic == False
        and clinvar_diff_aa_change.pathogenic == False
        and clinvar_same_nucleotide.pathogenic == False
        and clinvar_same_splice_site.pathogenic == False
    )


def test_splicing():
    """
    Test ClinVar annotation for splicing varaints
    """
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "test_var_splice_acceptor.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_files = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_files / file_name
    clinvar_dict = annotate_clinvar(
        variant.variant_info, variant.transcript_info, path_clinvar
    )
    for clinvar_type, clinvar_object in clinvar_dict.items():
        if clinvar_type.value == ClinVar_Type.SAME_AA_CHANGE.value:
            clinvar_same_aa_change = clinvar_object
        elif clinvar_type.value == ClinVar_Type.DIFF_AA_CHANGE.value:
            clinvar_diff_aa_change = clinvar_object
        elif clinvar_type.value == ClinVar_Type.SAME_NUCLEOTIDE.value:
            clinvar_same_nucleotide = clinvar_object
        elif clinvar_type.value == ClinVar_Type.SAME_SPLICE_SITE.value:
            clinvar_same_splice_site = clinvar_object
        else:
            continue
    assert (
        clinvar_same_aa_change.pathogenic == False
        and clinvar_diff_aa_change.pathogenic == False
        and clinvar_same_nucleotide.pathogenic == True
        and clinvar_same_nucleotide.ids == ["55502", "55501"]
        and clinvar_same_splice_site.pathogenic == True
        and clinvar_same_splice_site.ids == ["55502", "55501", "267593"]
    )
