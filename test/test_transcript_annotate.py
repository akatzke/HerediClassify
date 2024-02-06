#!/usr/bin/env python3

import pathlib

from test.example_variant_indel_GRCh38 import (
    create_example_del_frameshift,
    create_example_dup,
    create_example_del,
    create_example_ins,
)
from test.example_variant_splicing_GRCh38 import create_example_splice_acceptor_BRCA1
from test.example_variant_various_GRCh38 import create_start_lost
from test.example_variant_missense_GRCh38 import create_example_ter
from variant_classification.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from variant_classification.load_config import load_config
import test.paths as paths


def test_splicing():
    """
    Test annotation for a splice variant
    """
    test_trans, test_var = create_example_splice_acceptor_BRCA1()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_intronic.annotate(
        test_var,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        None,
        path_critical_region,
        test_trans,
    )
    assert (
        annot_trans.are_exons_skipped == True
        and annot_trans.is_NMD == True
        and annot_trans.is_truncated_region_disease_relevant == True
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_reading_frame_preserved == False
    )


def test_indel():
    """
    Test annotation for duplication
    """
    test_trans, test_var = create_example_dup()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        test_var,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        test_trans,
    )
    assert (
        annot_trans.is_NMD == True
        and annot_trans.is_reading_frame_preserved == False
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == True
        and round(annot_trans.diff_len_protein_percent, 2) == 0.29
        and annot_trans.frameshift == 1
        and annot_trans.ptc == 612
    )


def test_ter():
    """
    Test annotation for termination codon
    """
    test_trans, test_var = create_example_ter()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        test_var,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        test_trans[0],
    )
    assert (
        annot_trans.is_NMD == True
        and annot_trans.is_reading_frame_preserved == True
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == True
        and round(annot_trans.diff_len_protein_percent, 2) == 0.68
        and annot_trans.frameshift == 0
        and annot_trans.ptc == 275
    )


def test_del_inframe():
    """
    Test annotation for deletion
    """
    test_trans, test_var = create_example_del()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        test_var,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        test_trans,
    )
    assert (
        round(annot_trans.diff_len_protein_percent, 2) == 0.0
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == False
        and annot_trans.is_NMD == False
        and annot_trans.is_reading_frame_preserved == True
        and annot_trans.frameshift == 0
        and annot_trans.ptc == 392
    )


def test_del_frameshift():
    """
    Test annotation for deletion
    """
    test_trans, test_var = create_example_del_frameshift()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        test_var,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        test_trans,
    )
    assert (
        round(annot_trans.diff_len_protein_percent, 2) == -0.04
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == True
        and annot_trans.is_NMD == False
        and annot_trans.is_reading_frame_preserved == False
        and annot_trans.frameshift == -1
        and annot_trans.ptc == 565
    )


def test_ins():
    """
    Test annotation for insertion
    """
    test_trans, test_var = create_example_ins()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        test_var,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        test_trans,
    )
    assert (
        round(annot_trans.diff_len_protein_percent, 2) == 0.94
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == True
        and annot_trans.is_NMD == True
        and annot_trans.is_reading_frame_preserved == False
        and annot_trans.frameshift == 1
        and annot_trans.ptc == 113
    )


def test_start_loss():
    """
    Test annotation for a start loss variant
    """
    test_trans, test_var = create_start_lost()
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_region"]["root"]
    )
    path_critical_region = (
        dir_critical_region / config["annotation_files"]["critical_region"]["file"]
    )
    annot_trans = TranscriptInfo_start_loss.annotate(
        test_var, path_clinvar, path_uniprot, path_critical_region, test_trans
    )
    assert (
        round(annot_trans.diff_len_protein_percent, 2) == 0.05
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == True
        and annot_trans.exists_alternative_start_codon == True
        and annot_trans.position_alternative_start_codon
        == [35119569, 35119570, 35119571]
    )
