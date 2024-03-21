#!/usr/bin/env python3

import pathlib

from variant_classification.load_config import load_config
from variant_classification.check_disease_relevant_transcript import (
    check_disease_relevant_transcript,
)
from variant_classification.load_variant import load_variant
from variant_classification.ensembl import ensembl
from variant_classification.genotoscope_construct_variant_sequence import (
    construct_variant_coding_seq_exonic_variant,
)
from variant_classification.genotoscope_protein_len_diff import calculate_prot_len_diff
from variant_classification.genotoscope_assess_NMD import assess_NMD_threshold
from variant_classification.check_exon_disease_relevant import (
    check_exon_disease_relevant,
)

import test.paths as paths
from test.test_import_variant import create_json_string_from_variant


def test_NMD_threshold_true():
    """
    Test NMD threshold
    See that it is true
    """
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config = load_config(path_config)
    path_variant = paths.TEST / "test_variants" / "test_var_stop_gained_BRCA1_2.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
    var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
        transcript, variant_disease_relevant.variant_info, ref_transcript
    )
    diff_len_protein_percent, ptc = calculate_prot_len_diff(
        ref_transcript, var_seq, diff_len
    )
    NMD_threshold = {"ENST00000357654": 5418}
    try:
        nmd_threshold = NMD_threshold[transcript.transcript_id]
    except KeyError:
        raise KeyError("transcript not there")
    is_NMD, NMD_affected_exons = assess_NMD_threshold(
        transcript,
        variant_disease_relevant.variant_info,
        ptc,
        ref_transcript,
        diff_len,
        nmd_threshold,
    )
    assert is_NMD == True


def test_NMD_threshold_false():
    """
    Test NMD thresold
    See that it is False
    """
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config = load_config(path_config)
    path_variant = paths.TEST / "test_variants" / "test_var_stop_gained_BRCA1.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, config)
    transcript = variant_disease_relevant.transcript_info[0]
    ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
    var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
        transcript, variant_disease_relevant.variant_info, ref_transcript
    )
    NMD_threshold = {"ENST00000357654": 5418}
    diff_len_protein_percent, ptc = calculate_prot_len_diff(
        ref_transcript, var_seq, diff_len
    )
    try:
        nmd_threshold = NMD_threshold[transcript.transcript_id]
    except KeyError:
        raise KeyError("transcript not there")
    is_NMD, NMD_affected_exons = assess_NMD_threshold(
        transcript,
        variant_disease_relevant.variant_info,
        ptc,
        ref_transcript,
        diff_len,
        nmd_threshold,
    )
    assert is_NMD == False


def test_disease_relevant_transcript():
    """
    Test that filtering of disease relevant transcripts works
    """
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    final_config = load_config(path_config)
    path_variant = paths.TEST / "test_variants" / "test_var_stop_gained_BRCA1_2.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    diseases_relevant_transcripts = [
        entry["name"] for entry in final_config["disease_relevant_transcripts"]
    ]
    assert len(variant_disease_relevant.transcript_info) == len(
        diseases_relevant_transcripts
    ) and set(
        [transcript.transcript_id for transcript in variant.transcript_info]
    ) == set(
        diseases_relevant_transcripts
    )


def test_exon_not_disease_relevant():
    """
    Test that not disease relevant exon is recognised
    """
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    final_config = load_config(path_config)
    path_variant = paths.TEST / "test_variants" / "test_var_BRCA1_exon8.json"
    variant_str = create_json_string_from_variant(path_variant)
    variant = load_variant(variant_str)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    transcript = variant_disease_relevant.transcript_info[0]
    ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
    var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
        transcript, variant_disease_relevant.variant_info, ref_transcript
    )
    NMD_threshold = {"ENST00000357654": 5418}
    diff_len_protein_percent, ptc = calculate_prot_len_diff(
        ref_transcript, var_seq, diff_len
    )
    try:
        nmd_threshold = NMD_threshold[transcript.transcript_id]
    except KeyError:
        raise KeyError("transcript not there")
    is_NMD, NMD_affected_exons = assess_NMD_threshold(
        transcript,
        variant_disease_relevant.variant_info,
        ptc,
        ref_transcript,
        diff_len,
        nmd_threshold,
    )
    root_dir = pathlib.Path(final_config["annotation_files"]["root"])
    dir_critical_region = root_dir / pathlib.Path(
        final_config["annotation_files"]["critical_regions"]["root"]
    )
    path_disease_irrelevant_exons = (
        dir_critical_region
        / final_config["annotation_files"]["critical_regions"][
            "disease_irrelevant_exons"
        ]
    )
    is_affected_exon_disease_relevant = check_exon_disease_relevant(
        path_disease_irrelevant_exons, NMD_affected_exons
    )
    assert is_affected_exon_disease_relevant == False and is_NMD == True
