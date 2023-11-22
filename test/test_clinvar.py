#!/usr/bin/env python3

import pathlib

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
    path_config = pathlib.Path("./config.yaml")
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_files = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_files / file_name
    clinvar_dict = annotate_clinvar(test_var, [test_trans], path_clinvar)
    clinvar_same_aa_change = clinvar_dict[ClinVar_Type.SAME_AA_CHANGE]
    clinvar_diff_aa_change = clinvar_dict[ClinVar_Type.DIFF_AA_CHANGE]
    clinvar_same_nucleotide = clinvar_dict[ClinVar_Type.SAME_NUCLEOTIDE]
    clinvar_same_splice_site = clinvar_dict[ClinVar_Type.SAME_SPLICE_SITE]
    assert (
        clinvar_same_aa_change.pathogenic == True
        and clinvar_same_aa_change.ids == ["55147"]
        and clinvar_diff_aa_change.pathogenic == False
        and clinvar_same_nucleotide.pathogenic == False
        and clinvar_same_splice_site.pathogenic == False
    )


def test_splicing():
    """
    Test ClinVar annotation for splicing varaints
    """
    test_trans, test_var = create_example_splice_acceptor_BRCA1()
    path_config = pathlib.Path("./config.yaml")
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_files = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_files / file_name
    clinvar_dict = annotate_clinvar(test_var, [test_trans], path_clinvar)
    clinvar_same_aa_change = clinvar_dict[ClinVar_Type.SAME_AA_CHANGE]
    clinvar_diff_aa_change = clinvar_dict[ClinVar_Type.DIFF_AA_CHANGE]
    clinvar_same_nucleotide = clinvar_dict[ClinVar_Type.SAME_NUCLEOTIDE]
    clinvar_same_splice_site = clinvar_dict[ClinVar_Type.SAME_SPLICE_SITE]
    assert (
        clinvar_same_aa_change.pathogenic == False
        and clinvar_diff_aa_change.pathogenic == False
        and clinvar_same_nucleotide.pathogenic == True
        and clinvar_same_nucleotide.ids == ["55502", "55501", "55500"]
        and clinvar_same_splice_site.pathogenic == True
        and clinvar_same_splice_site.ids == ["55502", "55501", "55500"]
    )
