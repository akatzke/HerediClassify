#!/usr/bin/env python3

from test.test_import_variant import create_json_string_from_variant
import test.paths as paths

from variant_classification.classify import classify


def test_gene_specific_atm_frameshift():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "ATM_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_atm.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_atm_intron():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "ATM_intron_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_atm.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_atm_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "ATM_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_atm.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_atm_splice():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "ATM_splice_acceptor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_atm.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_atm_start_loss():
    """
    Test gene specific variant classification
    """
    path_variant = paths.TEST / "test_variants_gene_specific" / "ATM_start_lost.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_atm.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_frameshift():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA1_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_intron():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA1_intron_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_intron_2():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA1_intron_variant_2.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA1_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_splicing():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST
        / "test_variants_gene_specific"
        / "BRCA1_splice_acceptor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_stop_gained():
    """
    Test gene specific variant classification
    """
    path_variant = paths.TEST / "test_variants_gene_specific" / "BRCA1_stop_gained.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca1_synonymous():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA1_synonymous_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca2_frameshift():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA2_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca2_frameshift_2():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA2_frameshift_variant_2.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca2_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA2_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca2_splice():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "BRCA2_splice_donor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_brca2_stop_gained():
    """
    Test gene specific variant classification
    """
    path_variant = paths.TEST / "test_variants_gene_specific" / "BRCA2_stop_gained.json"
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_brca2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_cdh1_frameshift():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "CDH1_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_cdh1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_cdh1_frameshift_2():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "CDH1_frameshift_variant_2.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_cdh1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_cdh1_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "CDH1_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_cdh1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_cdh1_splice():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "CDH1_splice_donor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_cdh1.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_palb2_intron():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PALB2_intron_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_palb2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_palb2_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PALB2_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_palb2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_palb2_splice_acceptor():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST
        / "test_variants_gene_specific"
        / "PALB2_splice_acceptor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_palb2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_palb2_splice_donor():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PALB2_splice_donor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_palb2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_palb2_stop_gained():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PALB2_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_palb2.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_pten_frameshift():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PTEN_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_pten.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_pten_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PTEN_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_pten.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_pten_splice():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "PTEN_splice_donor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_pten.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_tp53_frameshift():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "TP53_frameshift_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_tp53.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_tp53_frameshift_2():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "TP53_frameshift_variant_2.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_tp53.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_tp53_missense():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "TP53_missense_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_tp53.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_tp53_missense_2():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "TP53_missense_variant_2.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_tp53.yaml"
    config_name, results = classify(path_config, var_str)


def test_gene_specific_tp53_splice():
    """
    Test gene specific variant classification
    """
    path_variant = (
        paths.TEST / "test_variants_gene_specific" / "TP53_splice_donor_variant.json"
    )
    var_str = create_json_string_from_variant(path_variant)
    path_config = paths.ROOT / "gene_specific" / "acmg_tp53.yaml"
    config_name, results = classify(path_config, var_str)
