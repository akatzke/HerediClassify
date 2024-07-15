#!/usr/bin/env

import pathlib
import logging
from functools import reduce, partial
import operator as op
from typing import Callable, Union, Any, Optional

from variant import Variant
from var_type import VARTYPE_GROUPS
import acmg_rules as Rules
from acmg_rules.utils import RuleResult, evidence_strength
from information import (
    Info,
    Classification_Info,
    Classification_Info_Groups,
)
from acmg_rules.computation_evidence_utils import (
    Threshold,
    THRESHOLD_DIRECTION,
)
from clinvar_annot import get_annotate_clinvar
from check_splice_site_classification_table import (
    get_annotate_splice_site_classification,
)
from check_splice_site_classification_table_include_last_exon_pos import (
    get_annotate_splice_site_classification_include_last_exon_pos,
)
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
    annotate_transcripts,
)
from check_coldspot_hotspot import (
    get_check_coldspot,
    get_check_hotspot,
)
from check_exon_pm5 import get_annotate_exon_classification_pm5
from check_splice_site_pm5_classification_table import (
    get_annotate_splice_site_classification_pm5,
)
from clinvar_annot_spliceai import (
    get_annotate_clinvar_spliceai_protein,
    get_annotate_clinvar_spliceai_splicing,
)
from clinvar_missense_similarity_score import (
    get_check_clinvar_missense_similarity,
)

logger = logging.getLogger("Classify.config_annotation")


def get_annotations_needed_from_rules(
    rule_list: list[str], class_info: Classification_Info
) -> dict[Callable, tuple[Info, ...]]:
    """
    Based on rule get Classification_Info objects required to apply the rules
    """
    RULE_DICTIONARY = {
        "pvs1": Rules.Pvs1,
        "pvs1_brca1": Rules.Pvs1_brca1,
        "pvs1_brca2": Rules.Pvs1_brca2,
        "pvs1_atm": Rules.Pvs1_atm,
        "pvs1_palb2": Rules.Pvs1_palb2,
        "pvs1_pten": Rules.Pvs1_pten,
        "pvs1_cdh1": Rules.Pvs1_cdh1,
        "ps1_protein": Rules.Ps1_protein,
        "ps1_protein_spliceai": Rules.Ps1_protein_spliceai,
        "ps1_protein_enigma": Rules.Ps1_protein_enigma,
        "ps1_splicing": Rules.Ps1_splicing,
        "ps1_splicing_clingen": Rules.Ps1_splicing_clingen,
        "ps1_protein_tp53": Rules.Ps1_protein_tp53,
        "ps1_splicing_tp53": Rules.Ps1_splicing_tp53,
        "ps1_splicing_pten": Rules.Ps1_splicing_pten,
        "ps3": Rules.Ps3,
        "pm1": Rules.Pm1,
        "pm1_supporting": Rules.Pm1_supporting,
        "pm1_tp53": Rules.Pm1_tp53,
        "pm2": Rules.Pm2,
        "pm2_supporting": Rules.Pm2_supporting,
        "pm2_supporting_faf": Rules.Pm2_supporting_faf,
        "pm2_supporting_less": Rules.Pm2_supporting_less,
        "pm2_supporting_less_faf": Rules.Pm2_supporting_less_faf,
        "pm2_supporting_no_indel": Rules.Pm2_supporting_no_ins_del_indel,
        "pm2_supporting_no_indel_faf": Rules.Pm2_supporting_no_ins_del_indel_faf,
        "pm4": Rules.Pm4,
        "pm4_pten": Rules.Pm4_pten,
        "pm4_stoploss": Rules.Pm4_stoploss,
        "pm5_protein": Rules.Pm5_protein,
        "pm5_protein_pathogenic": Rules.Pm5_protein_pathogenic,
        "pm5_protein_ptc": Rules.Pm5_protein_ptc,
        "pm5_splicing_ptc": Rules.Pm5_splicing_ptc,
        "pm5_protein_cdh1": Rules.Pm5_protein_cdh1,
        "pm5_splicing_cdh1": Rules.Pm5_splicing_cdh1,
        "pm5_enigma": Rules.Pm5_ptc_enigma,
        "pm5_protein_pten": Rules.Pm5_protein_pten,
        "pm5_protein_tp53": Rules.Pm5_protein_tp53,
        "pp1": Rules.Pp1,
        "pp2": Rules.Pp2,
        "pp3_splicing": Rules.Pp3_splicing,
        "pp3_splicing_enigma": Rules.Pp3_splicing_enigma,
        "pp3_splicing_mult_strength": Rules.Pp3_splicing_mult_strength,
        "pp3_splicing_cdh1": Rules.Pp3_splicing_cdh1,
        "pp3_protein": Rules.Pp3_protein,
        "pp3_protein_enigma": Rules.Pp3_protein_enigma,
        "pp3_protein_mult_strength": Rules.Pp3_protein_mult_strength,
        "pp4_enigma": Rules.Pp4_enigma,
        "ba1": Rules.Ba1,
        "ba1_faf": Rules.Ba1_faf,
        "ba1_with_absolute": Rules.Ba1_with_absolute,
        "bs1": Rules.Bs1,
        "bs1_faf": Rules.Bs1_faf,
        "bs1_with_absolute": Rules.Bs1_with_absolute,
        "bs1_supporting": Rules.Bs1_with_supporting,
        "bs1_supporting_faf": Rules.Bs1_with_supporting_faf,
        "bs1_absolute": Rules.Bs1_with_absolute,
        "bs2": Rules.Bs2,
        "bs2_supporting": Rules.Bs2_with_supporting,
        "bs3": Rules.Bs3,
        "bs4": Rules.Bs4,
        "bp1": Rules.Bp1,
        "bp1_annotation_cold_spot_strong": Rules.Bp1_annotation_cold_spot_strong,
        "bp3": Rules.Bp3,
        "bp4_splicing": Rules.Bp4_splicing,
        "bp4_splicing_enigma": Rules.Bp4_splicing_enigma,
        "bp4_splicing_mult_strength": Rules.Bp4_splicing_mult_strength,
        "bp4_protein": Rules.Bp4_protein,
        "bp4_protein_enigma": Rules.Bp4_protein_enigma,
        "bp4_protein_mult_strength": Rules.Bp4_protein_mult_strength,
        "bp5_enigma": Rules.Bp5_enigma,
        "bp7": Rules.Bp7,
        "bp7_deep_intronic_atm": Rules.Bp7_deep_intronic_atm,
        "bp7_deep_intronic_enigma": Rules.Bp7_deep_intronic_enigma,
        "bp7_deep_intronic_palb2": Rules.Bp7_deep_intronic_palb2,
    }
    rule_info_dict = {}
    for rule in rule_list:
        try:
            rule_class = RULE_DICTIONARY[rule.lower()]
        except KeyError:
            raise KeyError(
                f"{rule.lower()} not valid rule. \n Valid rules are {RULE_DICTIONARY.keys()}"
            )
        rule_fun, rule_args = rule_class.get_assess_rule(class_info)
        rule_info_dict[rule_fun] = rule_args
    return rule_info_dict


def get_unique_annotations_needed(
    fun_info_dict: dict[Callable, tuple[Info, ...]]
) -> list[Info]:
    """
    From dictionary mapping rules functions to their needed Classification_Info object
    """
    annots = []
    for rule_args in fun_info_dict.values():
        annots.append(rule_args)
    full_annot_list = [arg for args in fun_info_dict.values() for arg in args]
    unique_annot = [
        x for i, x in enumerate(full_annot_list) if i == full_annot_list.index(x)
    ]
    return unique_annot


def get_annotation_functions(
    annotations_needed: list[Info],
    variant: Variant,
    config: dict,
    class_info: Classification_Info,
) -> list[Info]:
    """
    Based on needed annotations perform annotation
    """
    ### Dictionary for all Classification_Info objects that are defined in the Variant object
    ### For definition Variant object see variant.py
    dict_annotation_variant = {
        class_info.VARIANT.name: lambda variant: partial(
            return_information, "variant_info", variant.variant_info
        ),
        class_info.TRANSCRIPT.name: lambda variant: partial(
            return_information, "transcript_info", variant.transcript_info
        ),
        class_info.VARIANT_CANCERHOTSPOTS.name: lambda variant: partial(
            return_information,
            "Cancer hotspots",
            variant.cancerhotspots,
        ),
        class_info.VARIANT_GNOMAD_POPMAX.name: lambda variant: partial(
            return_information, "GnomAD popmax", variant.gnomad_popmax
        ),
        class_info.VARIANT_GNOMAD_FAF.name: lambda variant: partial(
            return_information, "GnomAD faf", variant.gnomad_faf
        ),
        class_info.VARIANT_FLOSSIES.name: lambda variant: partial(
            return_information, "FLOSSIES", variant.flossies
        ),
        class_info.VARIANT_PREDICTION.name: lambda variant: partial(
            return_information, "Prediction tools", variant.prediction_tools
        ),
        class_info.FUNCTIONAL_ASSAY.name: lambda variant: partial(
            return_information, "Functional assay", variant.functional_assay
        ),
        class_info.SPLICING_ASSAY.name: lambda variant: partial(
            return_information, "Splicing assay", variant.splicing_assay
        ),
        class_info.VARIANT_MULTIFACTORIAL_LIKELIHOOD.name: lambda variant: partial(
            return_information,
            "Multifactorial likelihood",
            variant.multifactorial_likelihood,
        ),
    }

    ### Dictionary for all Classification_Info objects that have a get_annotation_function
    dict_annotation = {
        class_info.ANNOTATED_TRANSCRIPT_LIST.name: lambda variant, config: get_annotation_function_annotated_transcript(
            variant, config, class_info
        ),
        class_info.VARIANT_CLINVAR.name: lambda variant, config: get_annotation_function(
            get_annotate_clinvar, variant, config, class_info
        ),
        class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN.name: lambda variant, config: get_annotation_function(
            get_annotate_clinvar_spliceai_protein, variant, config, class_info
        ),
        class_info.VARIANT_CLINVAR_SPLICEAI_PROTEIN_SIMILARITY.name: lambda variant, config: get_annotation_function(
            get_check_clinvar_missense_similarity, variant, config, class_info
        ),
        class_info.VARIANT_CLINVAR_SPLICEAI_SPLICE.name: lambda variant, config: get_annotation_function(
            get_annotate_clinvar_spliceai_splicing, variant, config, class_info
        ),
        class_info.SPLICE_RESULT.name: lambda variant, config: get_annotation_function(
            get_annotate_splice_site_classification, variant, config, class_info
        ),
        class_info.SPLICE_RESULT_INCLUDE_LAST_EXON_POS.name: lambda variant, config: get_annotation_function(
            get_annotate_splice_site_classification_include_last_exon_pos,
            variant,
            config,
            class_info,
        ),
        class_info.VARIANT_HOTSPOT_ANNOTATION.name: lambda variant, config: get_annotation_function(
            get_check_hotspot, variant, config, class_info
        ),
        class_info.VARIANT_COLDSPOT_ANNOTATION.name: lambda variant, config: get_annotation_function(
            get_check_coldspot, variant, config, class_info
        ),
        class_info.SPLICE_RESULT_PM5.name: lambda variant, config: get_annotation_function(
            get_annotate_splice_site_classification_pm5, variant, config, class_info
        ),
        class_info.PM5_RESULTS_PTC.name: lambda variant, config: get_annotation_function(
            get_annotate_exon_classification_pm5, variant, config, class_info
        ),
    }

    ### Dictionary for all Classification_Info that belong to a Classification_Info_groups
    dict_annotation_groups = {
        Classification_Info_Groups.PATH: lambda annot, config: partial(
            get_path_from_config, annot.config_location, config
        ),
        Classification_Info_Groups.THRESHOLD_SINGLE: lambda annot, config: partial(
            get_threshold_from_config, annot.config_location, config
        ),
        Classification_Info_Groups.THRESHOLDS_LIKELIHOOD: lambda annot, config: partial(
            get_thresholds_likelihood,
            annot.config_location,
            config,
        ),
        Classification_Info_Groups.THRESHOLDS_PREDICTION: lambda annot, config: partial(
            get_thresholds_prediction,
            annot.config_location,
            config,
        ),
        Classification_Info_Groups.DISEASE_RELEVANT_TRANSCRIPT_THRESHOLD: lambda annot, config: partial(
            get_disease_relevant_transcript_thresholds, annot.config_location, config
        ),
        Classification_Info_Groups.CONFIG_ENTRY_STR: lambda annot, config: partial(
            get_config_entry_str, annot.config_location, config
        ),
    }

    for annotation in annotations_needed:
        if annotation.name in dict_annotation_variant.keys():
            annotation.compute_function = dict_annotation_variant[annotation.name](
                variant
            )
        elif annotation.name in dict_annotation.keys():
            annotation.compute_function = dict_annotation[annotation.name](
                variant, config
            )
        else:
            if annotation.group is Classification_Info_Groups.PATH:
                annotation.compute_function = dict_annotation_groups[
                    Classification_Info_Groups.PATH
                ](annotation, config)
            elif annotation.group is Classification_Info_Groups.THRESHOLD_SINGLE:
                annotation.compute_function = dict_annotation_groups[
                    Classification_Info_Groups.THRESHOLD_SINGLE
                ](annotation, config)
            elif annotation.group is Classification_Info_Groups.THRESHOLDS_LIKELIHOOD:
                annotation.compute_function = dict_annotation_groups[
                    Classification_Info_Groups.THRESHOLDS_LIKELIHOOD
                ](annotation, config)
            elif annotation.group is Classification_Info_Groups.THRESHOLDS_PREDICTION:
                annotation.compute_function = dict_annotation_groups[
                    Classification_Info_Groups.THRESHOLDS_PREDICTION
                ](annotation, config)
            elif (
                annotation.group
                is Classification_Info_Groups.DISEASE_RELEVANT_TRANSCRIPT_THRESHOLD
            ):
                annotation.compute_function = dict_annotation_groups[
                    Classification_Info_Groups.DISEASE_RELEVANT_TRANSCRIPT_THRESHOLD
                ](annotation, config)
            elif annotation.group is Classification_Info_Groups.CONFIG_ENTRY_STR:
                annotation.compute_function = dict_annotation_groups[
                    Classification_Info_Groups.CONFIG_ENTRY_STR
                ](annotation, config)
            else:
                raise ValueError(f"No annotation function defined for {annotation}.")
    return annotations_needed


def get_annotation_function_annotated_transcript(
    variant: Variant, config: dict, class_info: Classification_Info
) -> Callable[[], Any]:
    """
    Create annotation function for construction of Classification_Info.ANNOTATED_TRANSCIPT_LIST
    """
    relevant_classes = {
        VARTYPE_GROUPS.EXONIC: TranscriptInfo_exonic,
        VARTYPE_GROUPS.INTRONIC: TranscriptInfo_intronic,
        VARTYPE_GROUPS.START_LOST: TranscriptInfo_start_loss,
    }
    fun_dict = {}
    for name, entry in relevant_classes.items():
        fun_annot = prepare_function_for_annotation(
            partial(entry.get_annotate, class_info), variant, config, class_info
        )
        fun_dict[name] = fun_annot
    fun = partial(annotate_transcripts, variant, fun_dict)
    return fun


def get_annotation_function(
    get_fun: Callable, variant: Variant, config: dict, class_info: Classification_Info
):
    fun = prepare_function_for_annotation(
        partial(get_fun, class_info), variant, config, class_info
    )
    return fun


def prepare_function_for_annotation(
    fun: Callable[[], tuple[Callable, tuple[Info, ...]]],
    variant: Variant,
    config: dict,
    class_info: Classification_Info,
) -> Callable[[], Any]:
    """
    Prepare annotation function
    """
    annot_fun, annot_fun_args = fun()
    set_args = get_annotation_functions(
        list(annot_fun_args), variant, config, class_info
    )
    args = execute_annotation(set_args)
    for arg in args:
        if arg.value is None and not arg.optional:
            logger.warning(
                f"The annotation function {annot_fun} can not be defined, as {arg.name} is None. Annotation is skipped."
            )
            return lambda: None
    fun_annot = partial(annot_fun, *[argument.value for argument in args])
    return fun_annot


def execute_annotation(
    annotations_to_execute: list[Info],
) -> list[Info]:
    """
    Perform annotation for entrys in list
    If annotation does not execute, remove all rules that require that annotation from the fun_info_dict
    """
    for annotation in annotations_to_execute:
        if annotation.compute_function is not None:
            annotation.value = annotation.compute_function()
    return annotations_to_execute


def remove_rules_with_missing_annotation(
    fun_info_dict: dict[Callable, tuple[Info, ...]]
) -> dict[Callable, tuple[Info, ...]]:
    """
    Remove all rules that would try to access an information without a set value
    """
    original_dict = fun_info_dict.copy()
    for rule_fun, rule_args in original_dict.items():
        for annotation in rule_args:
            if annotation.value is None and not annotation.optional:
                del fun_info_dict[rule_fun]
                logger.info(f"Removed {rule_fun} from rules that will be assessed.")
                break
    return fun_info_dict


def apply_rules(fun_info_dict: dict[Callable, tuple[Info, ...]]) -> list[RuleResult]:
    """
    Execute all of the rules
    """
    rule_results = []
    for rule_fun, rule_args in fun_info_dict.items():
        rule_exec_fun = partial(rule_fun, *[arg.value for arg in rule_args])
        rule_result = rule_exec_fun()
        rule_results.append(rule_result)
    return rule_results


def return_information(info_name: str, info):
    """
    Return a given information
    """
    if info is None:
        logger.warning(f"{info_name} is not accessible.")
        return None
    return info


def get_path_from_config(
    config_location: Union[tuple[str, ...], None], config: dict
) -> Union[pathlib.Path, None]:
    """
    From config reconstruct full path
    """
    if config_location is None:
        logger.warning(
            f"No location in the configuration is defined for this path. Please check information.py."
        )
        return None
    try:
        root_files = pathlib.Path(config[config_location[0]]["root"])
        dir_files = root_files / pathlib.Path(
            config[config_location[0]][config_location[1]]["root"]
        )
        file_path = dir_files / pathlib.Path(
            config[config_location[0]][config_location[1]][config_location[2]]
        )
    except KeyError:
        logger.warning(
            f"The location {config_location} could not be found in the configuration file."
        )
        return None
    if not file_path.exists():
        logger.warning(
            f"The file {file_path} does not exist. Please make sure the file path is correct."
        )
        return None
    return file_path


def get_threshold_from_config(
    config_location: Optional[tuple[str, ...]], config: dict
) -> Optional[float]:
    """
    Get threshold from config
    """
    if config_location is None:
        logger.warning(
            f"No location in the configuration is defined for this threshold. Please check information.py."
        )
        return None
    try:
        config_value = reduce(op.getitem, config_location, config)
        try:
            assert isinstance(config_value, float) or isinstance(
                config_value, int
            ), f"Entry in configuration in location {' '.join(config_location)} is not of type float."
            return config_value
        except ValueError:
            logger.warning(
                f"The value at {config_location} is not a number. Please check your configuration."
            )
    except KeyError:
        logger.warning(
            f"The location {config_location} could not be found in the configuration file."
        )
    return None


def get_thresholds_likelihood(
    config_location: Union[tuple[str, ...], None], config: dict
) -> Union[Threshold, None]:
    """
    Get thresholds when multiple evidence strengths are defined
    """
    if config_location is None:
        logger.warning(
            f"No location in the configuration is defined for this threshold. Please check information.py."
        )
        return None
    try:
        config_prediction_tool: dict[str, Any] = reduce(
            op.getitem, config_location, config
        )
    except KeyError:
        logger.warning(
            f"The location {config_location} could not be found in configuration."
        )
        return None
    threshold_dict = dict(config_prediction_tool)
    try:
        dir = THRESHOLD_DIRECTION(threshold_dict["direction"])
    except ValueError:
        raise ValueError(
            f"The direction {threshold_dict['direction']} in {config_location} is not valid. Pleace check config."
        )
    del threshold_dict["direction"]
    thresholds = list(threshold_dict.values())
    assert all(
        isinstance(threshold, float) for threshold in thresholds
    ), f"Not all thresholds defined for {config_location[0]} are of type float. Please check."
    if "greater" in dir.value:
        keys_sorted = sorted(threshold_dict, key=threshold_dict.get, reverse=False)
    else:
        keys_sorted = sorted(threshold_dict, key=threshold_dict.get, reverse=True)
    strengths = [evidence_strength(key) for key in keys_sorted]
    return Threshold(
        name=config_location[0],
        direction=dir,
        thresholds=thresholds,
        strengths=strengths,
    )


def get_thresholds_prediction(
    config_location: Union[tuple[str, ...], None], config: dict
) -> Union[Threshold, None]:
    """
    Get thresholds when multiple evidence strengths are defined
    """
    if config_location is None:
        logger.warning(
            f"No location in the configuration is defined for this threshold. Please check information.py."
        )
        return None
    try:
        config_prediction_tool: dict[str, Any] = reduce(
            op.getitem, config_location, config
        )
    except KeyError:
        logger.warning(
            f"The location {config_location} could not be found in configuration."
        )
        return None
    threshold_dict = dict(config_prediction_tool)
    name_location = (*config_location[:-1], "name")
    name: str = reduce(op.getitem, name_location, config)
    try:
        dir = THRESHOLD_DIRECTION(threshold_dict["direction"])
    except ValueError:
        raise ValueError(
            f"The direction {threshold_dict['direction']} in {config_location} is not valid. Pleace check config."
        )
    del threshold_dict["direction"]
    thresholds = list(threshold_dict.values())
    if "greater" in dir.value:
        keys_sorted = sorted(threshold_dict, key=threshold_dict.get, reverse=False)
    else:
        keys_sorted = sorted(threshold_dict, key=threshold_dict.get, reverse=True)
    strengths = [evidence_strength(key) for key in keys_sorted]
    return Threshold(
        name=name,
        direction=dir,
        thresholds=thresholds,
        strengths=strengths,
    )


def get_disease_relevant_transcript_thresholds(
    config_location: Union[tuple[str, ...], None], config: dict
) -> Union[dict[str, int], None]:
    """
    Create dictionary containing the
    """
    if config_location is None:
        logger.warning(
            f"No location in the configuration is defined for this disease relevant transcript threshold. Please check information.py."
        )
        return None
    try:
        disease_relevant_transcripts = config["disease_relevant_transcripts"]
    except KeyError:
        logger.warning(
            f"No field disease_relevant_transcripts found in the configuration file."
        )
        return None
    try:
        thresh_dict = {}
        for transcript in disease_relevant_transcripts:
            thresh_dict[transcript["name"]] = transcript[config_location[-1]]
        if not bool(thresh_dict):
            return None
        return thresh_dict
    except KeyError:
        logger.warning(
            f"Either 'name' or {config_location[-1]} not defined for all disease relevant transcripts."
        )
        return None


def get_config_entry_str(
    config_location: Union[tuple[str, ...], None], config: dict
) -> Union[str, None]:
    """
    Create dictionary containing the
    """
    if config_location is None:
        logger.warning(f"No location in the configuration is defined.")
        return None
    try:
        conf_str = reduce(op.getitem, config_location, config)
        assert isinstance(
            conf_str, str
        ), f"Entry in configuration in location {' '.join(config_location)} is not of type string."
        return conf_str
    except KeyError:
        logger.warning(
            f"The location {config_location} could not be found in the configuration."
        )
        return None
