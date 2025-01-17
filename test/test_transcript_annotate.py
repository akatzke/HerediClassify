#!/usr/bin/env python3

import pathlib

from variant_classification.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)
from variant_classification.load_variant import load_variant
from variant_classification.check_disease_relevant_transcript import (
    check_disease_relevant_transcript,
)
from variant_classification.load_config import load_config
import test.paths as paths

from test.test_import_variant import create_json_string_from_variant


def test_splicing():
    """
    Test annotation for a splice variant
    """
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "test_var_splice_acceptor.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # Create relevant paths
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    # Execute transcript annotation
    annot_trans = TranscriptInfo_intronic.annotate(
        variant.variant_info,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        None,
        path_critical_region,
        transcript,
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
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "test_var_frameshift.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # Create relevant paths
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        variant.variant_info,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        transcript,
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
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "test_var_stop_gained_pms2.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # Create relevant paths
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    file_name_indel = "clinvar_small_indel.vcf.gz"
    path_clinvar_indel = dir_clinvar / file_name_indel
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        variant.variant_info,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        transcript,
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
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "test_var_inframe_deletion.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # Create relevant paths
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
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        variant.variant_info,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        transcript,
    )
    assert (
        round(annot_trans.diff_len_protein_percent, 2) == 0.0
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == False
        and annot_trans.is_NMD == False
        and annot_trans.is_reading_frame_preserved == True
        and annot_trans.frameshift == 0
        and annot_trans.ptc == 394
    )


def test_del_frameshift():
    """
    Test annotation for deletion
    """
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "CHEK2_stop_gained.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # test_trans, test_var = create_example_del_frameshift()
    # Create relevant paths
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
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        variant.variant_info,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        transcript,
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
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "BRCA1_frameshift_variant.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # test_trans, test_var = create_example_ins()
    # Create relevant paths
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
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        variant.variant_info,
        path_clinvar,
        path_clinvar_indel,
        path_uniprot,
        path_critical_region,
        None,
        None,
        transcript,
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
    # Load configuration
    path_config = paths.ROOT / "config.yaml"
    config = load_config(path_config)
    # Load transcript
    path_variant = paths.TEST / "test_variants" / "RAD51D_start_lost.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    # Create relevant paths
    root_dir = pathlib.Path(config["annotation_files"]["root"])
    dir_clinvar = root_dir / pathlib.Path(config["annotation_files"]["clinvar"]["root"])
    file_name = "clinvar_snv.vcf.gz"
    path_clinvar = dir_clinvar / file_name
    dir_uniprot = root_dir / pathlib.Path(config["annotation_files"]["uniprot"]["root"])
    path_uniprot = dir_uniprot / config["annotation_files"]["uniprot"]["rep"]
    dir_critical_region = root_dir / pathlib.Path(
        config["annotation_files"]["critical_regions"]["root"]
    )
    path_critical_region = (
        dir_critical_region
        / config["annotation_files"]["critical_regions"]["critical_region"]
    )
    annot_trans = TranscriptInfo_start_loss.annotate(
        variant.variant_info,
        path_clinvar,
        path_uniprot,
        path_critical_region,
        transcript,
    )
    assert (
        round(annot_trans.diff_len_protein_percent, 2) == 0.07
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_truncated_region_disease_relevant == True
        and annot_trans.exists_alternative_start_codon == True
        and annot_trans.position_alternative_start_codon
        == [35119569, 35119570, 35119571]
    )
