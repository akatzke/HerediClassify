#!/usr/bin/env python3


from test.test_import_variant import create_json_string_from_variant
import test.paths as paths

from variant_classification.ensembl import ensembl
from variant_classification.load_variant import load_variant
from variant_classification.load_config import load_config, get_gene_specific_config
from variant_classification.check_disease_relevant_transcript import (
    check_disease_relevant_transcript,
)
from variant_classification.genotoscope_construct_variant_sequence import (
    construct_variant_coding_seq_exonic_variant,
)
from variant_classification.genotoscope_protein_len_diff import (
    get_position_ptc,
    correct_position_ptc_for_indels,
)


def test_ptc_brca1():
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA1_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config = load_config(path_config)
    variant = load_variant(var_str)
    final_config = get_gene_specific_config(config, variant.variant_info.gene_name)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    transcript = variant_disease_relevant.transcript_info[0]
    ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
    var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
        transcript, variant_disease_relevant.variant_info, ref_transcript
    )
    codon_position_ptc = get_position_ptc(ref_transcript, var_seq.upper())
    corrected_codon_position_ptc = correct_position_ptc_for_indels(
        diff_len, codon_position_ptc
    )
    assert corrected_codon_position_ptc == 892


def test_ptc_brca2():
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA2_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config = load_config(path_config)
    variant = load_variant(var_str)
    final_config = get_gene_specific_config(config, variant.variant_info.gene_name)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    transcript = variant_disease_relevant.transcript_info[0]
    ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
    var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
        transcript, variant_disease_relevant.variant_info, ref_transcript
    )
    codon_position_ptc = get_position_ptc(ref_transcript, var_seq.upper())
    corrected_codon_position_ptc = correct_position_ptc_for_indels(
        diff_len, codon_position_ptc
    )
    assert corrected_codon_position_ptc == 24


def test_ptc_brca2_2():
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA2_frameshift_variant_2.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config = load_config(path_config)
    variant = load_variant(var_str)
    final_config = get_gene_specific_config(config, variant.variant_info.gene_name)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    transcript = variant_disease_relevant.transcript_info[0]
    ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
    var_seq, diff_len = construct_variant_coding_seq_exonic_variant(
        transcript, variant_disease_relevant.variant_info, ref_transcript
    )
    codon_position_ptc = get_position_ptc(ref_transcript, var_seq.upper())
    corrected_codon_position_ptc = correct_position_ptc_for_indels(
        diff_len, codon_position_ptc
    )
    assert corrected_codon_position_ptc == 459
