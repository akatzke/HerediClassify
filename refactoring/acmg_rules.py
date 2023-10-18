#!/usr/bin/env python3

from abc import ABC, abstractmethod
from typing import Callable

from refactoring.variant_annotate import Variant_annotated
from refactoring.transcript_annotated import *
from refactoring.rule_utils import (
    Prediction_result,
    aggregate_prediction_results,
    assess_prediction_tool,
    summarise_results_per_transcript,
)
from refactoring.config_annotation import classification_information


@dataclass
class RuleResult:
    name: str
    status: bool
    strength: str
    comment: str


class abstract_rule(ABC):

    @abstractmethod
    def get_assess_rule(cls) -> Callable:
        """
        Get function that assess rule
        """
        pass


    @abstractmethod
    def assess_rule(self, args) -> RuleResult:
        """
        Assess rule
        """
        pass

class pvs1(abstract_rule):
    """
    PVS1: Loss of function
    Devided into three separate parts: Frameshift, splice and start_loss
    """
    arguments = [classification_information.ANNOTATED_TRANSCRIPT_LIST]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule


    def assess_rule(self, annotated_transcripts: list[TranscriptInfo]) -> RuleResult:
        results = []
        for transcript in annotated_transcripts:
            if type(transcript) is TranscriptInfo_exonic:
                result_frameshift = self.assess_pvs1_frameshift_PTC(transcript)
                results.append(result_frameshift)
            elif type(transcript) is TranscriptInfo_intronic:
                result_splice = self.assess_pvs1_splice(transcript)
                results.append(result_splice)
            elif type(transcript) is TranscriptInfo_start_loss:
                result_start_loss = self.assess_pvs1_start_loss(transcript)
                results.append(result_start_loss)
        result = summarise_results_per_transcript(results)
        return result

    def assess_pvs1_start_loss(
        self, transcript: TranscriptInfo_start_loss
    ) -> RuleResult:
        """
        Assess PVS1 for start lost variants
        """
        if not transcript.exists_alternative_start_codon:
            comment = f"An alternative start code in transcript {transcript.transcript_id} at {transcript.position_alternative_start_codon} detected."
            result = RuleResult("PVS1_start_loss", False, "very_strong", comment)
        else:
            comment = f"No alternative start codon detected for transcript {transcript.transcript_id}."
            if transcript.is_truncated_region_disease_relevant:
                comment = f"Pathogenic variants detected between start codon and alternative start codon detected. \n ClinVar ID: {transcript.pathogenic_variants_truncated_region}"
                result = RuleResult("PVS1_start_loss", True, "moderate", comment)
            else:
                comment = "No pathogenic variant detected between start codon and alternative start codon."
                result = RuleResult("PVS1_start_loss", True, "supporting", comment)
        return result

    def assess_pvs1_splice(self, transcript: TranscriptInfo_intronic) -> RuleResult:
        """
        Assess PVS1 for splice variants
        """
        if transcript.is_NMD:
            comment = f"Transcript {transcript.transcript_id} undergoes NMD."
            if transcript.is_truncated_region_disease_relevant:
                comment = f"Skipped exon contais (likely) pathogenic variants and can therefore be considexoneredexon to be disease relevant. \n ClinVar ID: {transcript.pathogenic_variants_truncated_region}"
                result = RuleResult("PVS1_splice", True, "very_strong", comment)
            else:
                comment = "Skipped exon contains no (likely) pathogenic variants and is therefore not considered disease relevant."
                result = RuleResult("PVS1_splice", False, "very_strong", comment)
        elif (
            transcript.are_exons_skipped
            and not transcript.is_NMD
            and not transcript.is_reading_frame_preserved
        ):
            comment = f"Transcript {transcript.transcript_id} does not undergo NMD."
            if transcript.is_truncated_region_disease_relevant:
                comment = f"Skipped exon contais (likely) pathogenic variants and can therefore be considered to be disease relevant. \n ClinVar ID: {transcript.pathogenic_variants_truncated_region}"
                result = RuleResult("PVS1_splice", True, "strong", comment)
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Protein length change of"
                    result = RuleResult("PVS1_splice", True, "strong", comment)
                else:
                    comment = "Something"
                    result = RuleResult("PVS1_splice", True, "moderate", comment)
        elif transcript.are_exons_skipped and transcript.is_reading_frame_preserved:
            if transcript.is_truncated_region_disease_relevant:
                comment = "Something"
                result = RuleResult("PVS1_splice", True, "strong", comment)
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    result = RuleResult("PVS1_splice", True, "strong", comment)
                else:
                    comment = "Something"
                    result = RuleResult("PVS1_splice", True, "moderate", comment)
        else:
            comment = "Something"
            result = RuleResult("PVS1_splice", False, "very_strong", comment)
        return result

    def assess_pvs1_frameshift_PTC(
        self, transcript: TranscriptInfo_exonic
    ) -> RuleResult:
        """
        Assess PVS1 for frameshift variants
        """
        if transcript.is_NMD:
            if transcript.is_truncated_region_disease_relevant:
                comment = "Something"
                result = RuleResult("PVS1_splice", True, "very_strong", comment)
            else:
                comment = "Something"
                result = RuleResult("PVS1_splice", False, "very_strong", comment)
        else:
            if transcript.is_truncated_region_disease_relevant:
                comment = "Something"
                result = RuleResult("PVS1_splice", True, "strong", comment)
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    result = RuleResult("PVS1_splice", True, "strong", comment)
                else:
                    comment = "Something"
                    result = RuleResult("PVS1_splice", True, "moderate", comment)
        return result


class ps1(abstract_rule):
    """
    PS1: Same amino accid change has been identified as pathogenic in ClinVar
    """
    arguments = [classification_information.VARIANT_CLINVAR]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(self, variant: Variant_annotated) -> RuleResult:
        if variant.same_aa_change_clinvar.pathogenic:
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {variant.same_aa_change_clinvar.ids}."
            result = RuleResult("PS1", True, "strong", comment)
        else:
            comment = "No matches found for variant."
            result = RuleResult("PS1", False, "strong", comment)
        return result


class pm1(abstract_rule):
    """
    PM1: Variant located in mutational hot spot or citical protein region
    """
    arguments = [classification_information.VARIANT_ANNOT]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(self, variant: Variant_annotated) -> RuleResult:
        if variant.affected_region.critical_region:
            comment = f"Variant in mutational hotspot. {variant.affected_region.critical_region_type}"
            result = RuleResult("PM1", True, "moderate", comment)
        else:
            comment = "Variant not in mutational hotspot"
            result = RuleResult("PM1", False, "moderate", comment)
        return result


class pm2(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    """
    arguments = [classification_information.VARIANT_ANNOT, classification_information.THRESHOLD_PM2]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule


    def assess_rule(self, variant: Variant_annotated, threshold_pm2: float) -> RuleResult:
        if variant.gnomad.frequency > threshold_pm2:
            comment = f"Variant occures with {variant.gnomad.frequency} in {variant.gnomad.name}."
            result = RuleResult("PM2", False, "moderate", comment)
        else:
            comment = f"Variant occurs with {variant.gnomad.frequency} in {variant.gnomad.name}."
            result = RuleResult("PM2", True, "moderate", comment)
        return result


class pm4(abstract_rule):
    """
    PM4: Protein length change caused by variant is above 10% threshold
    """
    arguments = [classification_information.ANNOTATED_TRANSCRIPT_LIST, classification_information.THRESHOLD_DIFF_LEN_PROT_PERCENT]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(
            self, annotated_transcript_list: list[TranscriptInfo_annot], threshold_diff_len_prot_percent: float
    ) -> RuleResult:
        results = []
        for transcript in annotated_transcript_list:
            if not transcript.is_truncated_region_disease_relevant:
                comment = (
                    f"Transcript {transcript.transcript_id} is not disease relevant."
                )
                result = RuleResult("PM4", False, "moderate", comment)
            elif (
                transcript.diff_len_protein_percent
                > threshold_diff_len_prot_percent
                and transcript.is_truncated_region_disease_relevant
                and not transcript.len_change_in_repetitive_region
            ):
                comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region does not overlap repetitive region."
                result = RuleResult("PM4", True, "moderate", comment)
            elif (
                transcript.diff_len_protein_percent > 0.1
                and transcript.is_truncated_region_disease_relevant
                and transcript.len_change_in_repetitive_region
            ):
                comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
                result = RuleResult("PM4", False, "moderate", comment)
            else:
                comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}"
                result = RuleResult("PM4", False, "moderate", comment)
            results.append(result)
        final_result = summarise_results_per_transcript(results)
        return final_result


class pm5(abstract_rule):
    """
    PM5: Pathogenic missense variant to different amino acid in same position classified as pathogenic in ClinVar
    """
    arguments = [classification_information.VARIANT_CLINVAR]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(self, variant: Variant_annotated) -> RuleResult:
        if variant.diff_aa_change_clinvar.pathogenic:
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {variant.diff_aa_change_clinvar.pathogenic}."
            result = RuleResult("PM5", True, "moderate", comment)
        else:
            comment = "No matches found for variant."
            result = RuleResult("PM5", False, "moderate", comment)
        return result


class bp4_pp3(abstract_rule):
    """
    BP4 and PP3: Assess results of prediction programs
    """
    arguments = [classification_information.VARIANT_ANNOT, classification_information.THRESHOLD_PATHOGENICITY_PREDICTION, classification_information.THRESHOLD_SPLICING_PREDICTION]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(self, variant: Variant_annotated) -> RuleResult:
        splicing_prediction = []
        pathogenicity_prediction = []
        conservation_prediction = []
        for entry in variant.prediction_tools:
            if entry.type == "pathogenicity":
                result = assess_prediction_tool(entry, variant.thresholds[entry.name])
                pathogenicity_prediction.append(result)
            elif entry.type == "splcing":
                result = assess_prediction_tool(entry, variant.thresholds[entry.name])
                splicing_prediction.append(result)
            elif entry.type == "conservation":
                result = assess_prediction_tool(entry, variant.thresholds[entry.name])
                conservation_prediction.append(result)
            else:
                print("There is an error here")
                break
        splicing_result = aggregate_prediction_results(splicing_prediction)
        pathogenicity_result = aggregate_prediction_results(pathogenicity_prediction)
        is_conserved = all(conservation_prediction)
        if (
            splicing_result == Prediction_result.PATHOGENIC
            and pathogenicity_result == Prediction_result.PATHOGENIC
            and is_conserved
        ):
            comment = "Varinat is predicted to affect splicing and to be pathogenic and is conserved."
            result = RuleResult("PP3", True, "supporting", comment)
            result = RuleResult("BP4", False, "supporting", comment)
        elif (
            splicing_result == Prediction_result.BENIGN
            and pathogenicity_result == Prediction_result.BENIGN
            and not is_conserved
        ):
            comment = "Varinat is predicted to not affect splicing and to not be pathogenic and is not conserved."
            result = RuleResult("PP3", False, "supporting", comment)
            result = RuleResult("BP4", True, "supporting", comment)
        elif splicing_result == Prediction_result.PATHOGENIC and is_conserved:
            comment = "Variant is predicted to affect splicing."
            result = RuleResult("PP3", True, "supporting", comment)
            result = RuleResult("BP4", False, "supporting", comment)
        elif pathogenicity_result == Prediction_result.PATHOGENIC and is_conserved:
            comment = "Variant is predicted to be pathogenic."
            result = RuleResult("PP3", True, "supporting", comment)
            result = RuleResult("BP4", False, "supporting", comment)
        elif splicing_result == Prediction_result.BENIGN and not is_conserved:
            comment = "Variant is predicted to not affect splicing and is not coserved."
            result = RuleResult("PP3", False, "supporting", comment)
            result = RuleResult("BP4", True, "supporting", comment)
        elif pathogenicity_result == Prediction_result.BENIGN and not is_conserved:
            comment = "Variant is predicted to be pathogenic."
            result = RuleResult("PP3", False, "supporting", comment)
            result = RuleResult("BP4", True, "supporting", comment)
        else:
            comment = "No conclusive evidence regarding the prediction of pathogenicity and splicing."
            result = RuleResult("PP3", False, "supporting", comment)
            result = RuleResult("BP4", False, "supporting", comment)
        return result


class ba1(abstract_rule):
    """
    BA1: High frequency of variant in healthy population (e.g. gnomAD)
    """
    arguments = [classification_information.VARIANT_ANNOT, classification_information.THRESHOLD_BA1]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(self, variant: Variant_annotated, threshold_ba1: float) -> RuleResult:
        if variant.gnomad.frequency > threshold_ba1:
            comment = f"Variant occures with {variant.gnomad.frequency} in {variant.gnomad.name}."
            result = RuleResult("BA1", False, "stand_alone", comment)
        else:
            comment = f"Variant occurs with {variant.gnomad.frequency} in {variant.gnomad.name}."
            result = RuleResult("BA1", True, "stand_alone", comment)
        return result


class bs1(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """
    arguments = [classification_information.VARIANT_ANNOT, classification_information.THRESHOLD_BS1]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_bs1(self, variant: Variant_annotated, threshold_bs1: float) -> RuleResult:
        if variant.gnomad.frequency > threshold_bs1:
            comment = f"Variant occures with {variant.gnomad.frequency} in {variant.gnomad.name}."
            result = RuleResult("BS1", False, "strong", comment)
        else:
            comment = f"Variant occurs with {variant.gnomad.frequency} in {variant.gnomad.name}."
            result = RuleResult("BS1", True, "strong", comment)
        return result


class bs2(abstract_rule):
    """
    BS2: Mutation found in healthy individuals
    """
    arguments = [classification_information.VARIANT_ANNOT, classification_information.THRESHOLD_BS2]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_rule(self, variant: Variant_annotated, threshold_bs2: float) -> RuleResult:
        if variant.flossies.frequency == threshold_bs2:
            comment = "Something"
            result = RuleResult("BS2", True, "strong", comment)
        else:
            comment = "Something"
            result = RuleResult("BS2", False, "strong", comment)
        return result


class bp3(abstract_rule):
    """
    BP3: Protein length change in repetitive region
    """
    arguments = [classification_information.ANNOTATED_TRANSCRIPT_LIST, classification_information.THRESHOLD_DIFF_LEN_PROT_PERCENT]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_bp3(self, annotated_transcripts: list[TranscriptInfo], threshold_diff_len_prot_percent: float) -> RuleResult:
        results = []
        for transcript in annotated_transcripts:
            if (
                type(transcript) != TranscriptInfo_exonic
                or type(transcript) != TranscriptInfo_intronic
            ):
                break
            if not transcript.transcript_disease_relevant:
                comment = (
                    f"Transcript {transcript.transcript_id} is not disease relevant."
                )
                result = RuleResult("BP3", False, "supporting", comment)
            elif (
                transcript.diff_len_protein_percent
                <= threshold_diff_len_prot_percent
                and transcript.len_change_in_repetitive_region
            ):
                comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
                result = RuleResult("BP3", True, "supporting", comment)
            else:
                comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}."
                result = RuleResult("BP3", False, "supporting", comment)
            results.append(result)
        final_result = summarise_results_per_transcript(results)
        return final_result


class bp7(abstract_rule):
    """
    BP7: Silent missense variant is predicted to have effect on splicing
    """
    arguments = [classification_information.VARIANT_ANNOT, classification_information.THRESHOLD_SPLICING_PREDICTION]

    @classmethod
    def get_assess_rule(cls) -> Callable:
        return cls.assess_rule

    def assess_bp7(self, variant: Variant_annotated) -> RuleResult:
        splicing_prediction = []
        conservation_prediction = []
        for entry in variant.prediction_tools:
            if entry.type == "splcing":
                result = assess_prediction_tool(entry, variant.thresholds[entry.name])
                splicing_prediction.append(result)
            elif entry.type == "conservation":
                result = assess_prediction_tool(entry, variant.thresholds[entry.name])
                conservation_prediction.append(result)
            else:
                print("There is an error here")
                break
        splicing_result = aggregate_prediction_results(splicing_prediction)
        is_conserved = all(conservation_prediction)
        if splicing_result == Prediction_result.BENIGN and not is_conserved:
            comment = (
                "Varinat is predicted to not affect splicing and is not conserved."
            )
            result = RuleResult("BP4", True, "supporting", comment)
        elif splicing_result == Prediction_result.BENIGN and is_conserved:
            comment = "Variant is predicted to not affect splicing, but is coserved."
            result = RuleResult("BP4", False, "supporting", comment)
        else:
            comment = "No conclusive evidence regarding the prediction of splicing."
            result = RuleResult("BP4", False, "supporting", comment)
        return result
