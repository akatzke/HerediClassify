#!/usr/bin/env python3

import pathlib

import test.paths as paths
from test.example_variant_missense_GRCh38 import create_example_mis_BRCA1_2
from test.example_variant_splicing_GRCh38 import create_example_splice_acceptor_BRCA1
from variant_classification.clinvar_annot import annotate_clinvar
from variant_classification.clinvar_utils import ClinVar_Type
from variant_classification.load_config import load_config


def test_missense_2():
    """
    Test ClinVar annotation for missense variant
    """
    test_trans, test_var = create_example_mis_BRCA1_2()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_files = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_files / file_name
    clinvar_dict = annotate_clinvar(test_var, [test_trans], path_clinvar)
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
    test_trans, test_var = create_example_splice_acceptor_BRCA1()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_files = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_files / file_name
    clinvar_dict = annotate_clinvar(test_var, [test_trans], path_clinvar)
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
        and clinvar_same_splice_site.ids == ["55502", "55501"]
    )
