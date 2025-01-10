Implemented rules
^^^^^^^^^^^^^^^^^^^^^^^^

The decision trees for all implemented criteria are listed below.
Please also check the gene specific recommendations for more details.

PVS1
~~~~
Each implementation of PVS1 implements a decision tree to determine whether or not a variant causes a loss of function for the relevant protein.

For further information for the general implementation of PVS1 see Abou Tayoun et al. (Human Mutation, 2018).
The gene specific implementations for PVS1 follow decision trees based on that suggested by Abou Tayoun, and specified in the corresponding gene-specific guidelines on ClinGen.

- **PVS1**
- **PVS1_BRCA1**
- **PVS1_BRCA2**
- **PVS1_ATM**
- **PVS1_PALB2**
- **PVS1_PTEN**
- **PVS1_CDH1**

PS1
~~~~
Check ClinVar for variants with the same amino acid exchange and whether or not they are (likely) pathogenic.
For specifications of the individual rules please see list of implementations below.


- **PS1_protein**:
  Check for variant with same amino acid exchange in ClinVar.
  ClinVar classification has to be pathogenic or likely pathogenic/pathogenic in order for a variant to be considered for PS1 evidence.
  A likely pathogenic classification on ClinVar is not sufficient.

- **PS1_protein_spliceai**:
  Same as PS1_protein.
  In addition, the variant under assessment and ClinVar variants must be predicted to not affect splicing.

- **PS1_protein_enigma**:
  Same as PS1_protein_spliceai.
  In case of only likely pathogenic variants with the same amino acid exchange, the evidence strength of PS1 is downgraded to moderate.

- **PS1_splicing**:
  Variant at same position has been classified as pathogenic.

- **PS1_splicing_clingen**:
  Implementing the decision tree for splicing variants, as recommended by Walker et al (AJHG, 2023).

- **PS1_protein_tp53**:
  See recommendation for the application of PS1 for TP53 on ClinGen.

- **PS1_splicing_tp53**:
  See recommendation for the application of PS1 for TP53 on ClinGen.

- **PS1_splicing_pten**:
  See recommendations for the application of PS1 for PTEN on ClinGen

PS3
~~~
Functional assay indicated pathogenicity of variant.

- PS3

PM1
~~~
Variant located in mutational hotspot.

- **PM1_supporting**:
  PM1 is only implemented based on supporting evidence strength.
  For the identification of hotspot regions, see download scripts.
  In the configuration file under annotation_files/critical_region/critical_regions a bed file with defined critical regions can be defined that is being used for PM1.

- **PM1_TP53**:
  See recommendations for the application of PM1 for TP53


PM2
~~~
Variant absent from healthy population.

- **PM2**:
  Applies PM2 with moderate evidence strength.
  PM2 will be using the popmax allele frequency and fall back on general allele frequency should popmax not be available.
  Allele frequency has to be less than or equal to threshold set for PM2.

- **PM2_supporting_faf**:
  Uses filtered allele frequency instead of popmax allele frequeny.

- **PM2_supporting_less**:
  When comparing allele frequency to allele frequency, the allele frequency has to be smaller than the threshold.

- **PM2_supporting_less_faf**:
  When comparing allele frequency to allele frequency, the allele frequency has to be smaller than the threshold.
  Uses filtered allele frequency instead of popmax allele frequeny.

- **PM2_supporting_no_indel**:
  PM2 will not be applied to indels.

- **PM2_supporting_no_indel_faf**:
  PM2 will not be applied to indels.
  Uses filtered allele frequency instead of popmax allele frequeny.

PM4
~~~
Length change of protein larger than 10% and located in disease relevant region.

- **PM4**:
  PM4 applies to inframe insertions or deletions, stop loss variants, and frameshift variants.
  For inframe variants, deletion/insertion has to be either located in a disease relevant region, or the change in protein length must be greater than 10%.
  PM4 is applied to all stop loss variants that extend the protein by more than 10%. In case a repetitive region is affected PM4 does not apply.
  PM4 applies to frameshift variants in case protein length is changed by more than 10% and no repetitive region is affected.

- **PM4_stoploss**:
  PM4 only applies to stop loss variants that increase protein length by more than 10% and affect a disease relevant region.

- **PM4_PTEN**:
  See recommendations for the application of PM4 for PTEN.
  PM4 applies to inframe deletions and insertions if they are located in disease relevant region and to all increases in protein length if not located in repetitive regions.


PM5
~~~
ClinVar entry for same amino acid position with different amino acid exchange has been shown to be pathogenic.

- **PM5_protein**:
  PM5 applies if a variant with different amino acid exchange at same position has been classifed as pathogenic, pathogenic/likely pathogenic or likely pathogenic in ClinVar.

- **PM5_protein_pathogenic**:
  PM5 applies if a variant with different amino acid exchange at same position has been classifed as pathogenic or likely pathogenic/ pathogenic in ClinVar.
  A likely pathogenic classification on ClinVar is not sufficient.

- **PM5_protein_ptc**:
  PM5 is applicable if the premature termination codon (PTC) caused by the variant is located upstream of the last known pathogenic PTC.
  This applies to variant types: frameshift, stop gained and stop loss.
  The last known PTC has to be provided in the configuration file under disease_relevant_transcripts under pos_last_known_patho_ptc (int). Please also specify which transcript (Ensemble) the PTC is located in.
  The last kown PTC can be defined for multiple transcripts.

- **PM5_splicing_ptc**:
  PM5 is applicable if the premature termination codon (PTC) caused by the variant is located upstream of the last known pathogenic PTC.
  This applies to variant types: splice acceptor variant and splice donor variant.
  The last known PTC has to be provided in the configuration file under disease_relevant_transcripts under pos_last_known_patho_ptc (int). Please also specify which transcript (Ensemble) the PTC is located in.
  The last kown PTC can be defined for multiple transcripts.

- **PM5_protein_CDH1**:
  PM5 is applicable if the premature termination codon (PTC) caused by the variant is located upstream of the last known pathogenic PTC.
  Alternatively, PM5 applies when a variant is predicted to cause nonsense mediated decay.
  This applies to variant types: frameshift, stop gained and stop loss.
  The last known PTC has to be provided in the configuration file under disease_relevant_transcripts under pos_last_known_patho_ptc (int). Please also specify which transcript (Ensemble) the PTC is located in.

- **PM5_splicing_CDH1**:
  PM5 is applicable if the premature termination codon (PTC) caused by the variant is located upstream of the last known pathogenic PTC.
  Alternatively, PM5 applies when a variant is predicted to cause nonsense mediated decay.
  This applies to variant types: splice acceptor variant and splice donor variant.
  The last known PTC has to be provided in the configuration file under disease_relevant_transcripts under pos_last_known_patho_ptc (int). Please also specify which transcript (Ensemble) the PTC is located in.

- **PM5_enigma**:
  PM5 applies based on the exon the premature termination codon is located in.
  A table with the necessary data is provided under data/PM5_annotations/PM5_PTC_BRCA1.csv

- **PM5_protein_PTEN**:
  PM5 applies if a variant with different amino acid exchange at same position has been classified as pathogenic, pathogenic/likely pathogenic or likely pathogenic in ClinVar.
  A splicing effect has to be excluded for the variant under assessment as well as for the variants in ClinVar.
  Additionally, the blosum62 score of the variant under assessments needs to be smaller than that of the ClinVar variants.

- **PM5_protein_TP53**:
  PM5 applies if a variant with different amino acid exchange at same position has been classified as pathogenic, pathogenic/likely pathogenic or likely pathogenic in ClinVar.
  A splicing effect has to be excluded for the variant under assessment as well as for the variants in ClinVar.
  Additionally, the Grantham score of the variant under assessments needs to be greater than that of the ClinVar variants.

PP1
~~~
Variant segregates with disease.

- **PP1**:
  Variant segregates with disease.
  Threshold for likelihood can be set under likelihood_thresholds for benignity and pathogenicity and differnet evidence strengths.

PP2
~~~
Gene is known to have any pathogenic variants of the same type.

- **PP2**:
  PP2 applies to all missense variants.

PP3
~~~
Computational evidence for pathogenicity of variant.

- **PP3_splicing**:
  Checks if variant is predicted to be pathogenic by prediction tool.
  Threshold can be set under prediction_tool_thresholds/splicing_prediction/pathogenic/supporting (type:float).

- **PP3_splicing_enigma**:
  Checks if variant is predicted to be pathogenic by prediction tool.
  Threshold can be set under prediction_tool_thresholds/splicing_prediction/pathogenic/supporting (type:float).
  PP3 can not be applied if variant is located outside of disease relevant regions.

- **PP3_splicing_enigma_mult_strength**:
  PP3 can not be applied if variant is located outside of disease relevant regions.
  Thresholds to be used for PP3 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/splicing_prediction/pathogenic/(supporting or moderate or strong or very_strong) (type:float).

- **PP3_splicing_mult_strength**:
  Thresholds to be used for PP3 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/splicing_prediction/pathogenic/(supporting or moderate or strong or very_strong) (type:float).

- **PP3_splicing_CDH1**:
  PP3 only applies to intronic variants located outside of the canonical splice site.
  Checks if variant is predicted to be pathogenic by prediction tool.
  Threshold can be set under prediction_tool_thresholds/splicing_prediction/pathogenic/supporting (type:float).

- **PP3_protein**:
  Checks if variant is predicted to be pathogenic by prediction tool.
  Threshold can be set under prediction_tool_thresholds/pathogenicity_prediction/pathogenic/supporting (type:float).

- **PP3_protein_enigma**:
  Checks if variant is predicted to be pathogenic by prediction tool.
  Threshold can be set under prediction_tool_thresholds/pathogenicity_prediction/pathogenic/supporting (type:float).
  PP3 can not be applied if variant is located outside of disease relevant regions.

- **PP3_protein_enigma_mult_strength**:
  Thresholds to be used for PP3 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/pathogenicity_prediction/pathogenic/(supporting or moderate or strong or very_strong) (type:float).
  PP3 can not be applied if variant is located outside of disease relevant regions.

- **PP3_protein_mult_strength**:
  Thresholds to be used for PP3 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/pathogenicity_prediction/pathogenic/(supporting or moderate or strong or very_strong) (type:float).


BA1
~~~
Variant is very common in helathy popualtion.

- **BA1**:
  BA1 is using the popmax allele frequency and falls back on genereal allele frequency should popmax not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_ba1.

- **BA1_faf**:
  BA1 is using the filtered allele frequency and falls back on genereal allele frequency should filtered allele frequency not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_ba1.

- **BA1_with_absolute**:
  BA1 is using the popmax allele frequency and falls back on genereal allele frequency should popmax not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_ba1.
  Additionally, absolute allele count in popmax allele count is checked.
  Threshold can be set under allele_frequency_thresholds/threshold_ba1_absolute.

BS1
~~~
Variant is common in helathy popualtion.

- **BS1**:
  BS1 is using the popmax allele frequency and falls back on general allele frequency should popmax not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_bs1 (type:float).

- **BS1_faf**:
  BS1 is using the filtered allele frequency and falls back on general allele frequency should filtered allele frequency not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_bs1 (type:float).

- **BS1_with_absolute**:
  BS1 is using the popmax allele frequency and falls back on general allele frequency should popmax not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_bs1 (type:float).
  Additionally, absolute allele count in popmax allele count is checked.
  Threshold can be set under allele_frequency_thresholds/threshold_bs1_absolute (type:int).

- **BS1_supporting**:
  BS1 is using the popmax allele frequency and falls back on general allele frequency should popmax not be available.
  Threshold can be set under allele_frequency_thresholds/threshold_bs1 (type:float) and allele_frequency_thresholds/threshold_bs1_supporting (type:float).

- **BS1_supporting_faf**:
  BS1 is using the filtered allele frequency and falls back on general allele frequency should filtered allele frequency not be available.
  Option to apply BS1 with strong and supporting evidence strength.
  Threshold can be set under allele_frequency_thresholds/threshold_bs1 (type:float) and allele_frequency_thresholds/threshold_bs1_supporting (type:float).
  Both thresholds need to be given.

BS2
~~~
Mutation found in healthy individual.

- **BS2**:
  Checks FLOSSIES database for presence of this variant.
  Threshold can be set under allele_frequency_thresholds/threshold_bs2 (type:int).

- **BS2_supporting**:
  Option to apply BS2 with strong and supporting evidence strength.
  Threshold can be set under allele_frequency_thresholds/threshold_bs2 (type:int) and allele_frequency_thresholds/threshold_bs2_supporting.
  Both thresholds need to be given.

BS3
~~~
Functional data indicating benignity.

- **BS3**


BS4
~~~
Variant does not segregate with disease.

- **BS4**:
  Variant does not segregate with disease.
  Threshold for likelihood can be set under likelihood_thresholds for benignity and pathogenicity and differnet evidence strengths.

BP1
~~~
Missense variant in a gene where missense variants are known not to be disease causative.

- **BP1**:
  BP1 applies to all missense variants.
- **BP1_annotation_cold_spot_strong**:
  Variant located in known cold spot region.
  BP1 is applied with strong evidence strength.
  Bed file with cold-spot regions can be defined under annotation_files/critical_regions/coldspot_region.


BP3
~~~
Variant located in repetitive region.

- **BP3**:
  Check if variant causes differen in portein length and if length change is located in repetitive region, BP3 is applied.

BP4
~~~
Computational evidence for benignity of variant.

- **BP4_splicing**:
  Checks if variant is predicted to be benign by prediction tool.
  Threshold can be set under prediction_tool_thresholds/splicing_prediction/benign/supporting (type:float).

- **BP4_splicing_enigma**:
  Checks if variant is predicted to be benign by prediction tool.
  Threshold can be set under prediction_tool_thresholds/splicing_prediction/benign/supporting (type:float).
  BP4 can not be applied if variant is located outside of disease relevant regions.

- **BP4_splicing_enigma_mult_strength**:
  BP4 can not be applied if variant is located outside of disease relevant regions.
  Thresholds to be used for BP4 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/splicing_prediction/benign/(supporting or moderate or strong or very_strong) (type:float).

- **BP4_splicing_mult_strength**:
  Thresholds to be used for BP4 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/splicing_prediction/benign/(supporting or moderate or strong or very_strong) (type:float).

- **BP4_protein**:
  Checks if variant is predicted to be benign by prediction tool.
  Threshold can be set under prediction_tool_thresholds/pathogenicity_prediction/benign/supporting (type:float).

- **BP4_protein_enimga**:
  Checks if variant is predicted to be benign by prediction tool.
  Threshold can be set under prediction_tool_thresholds/pathogenicity_prediction/benign/supporting (type:float).
  BP4 can not be applied if variant is located outside of disease relevant regions.

- **BP4_protein_enimga_mult_strength**:
  Thresholds to be used for BP4 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/pathogenicity_prediction/benign/(supporting or moderate or strong or very_strong) (type:float).
  BP4 can not be applied if variant is located outside of disease relevant regions.

- **BP4_protein_mult_strength**:
  Thresholds to be used for BP4 for all evidence strength levels are assessed and can be set under prediction_tool_thresholds/pathogenicity_prediction/benign/(supporting or moderate or strong or very_strong) (type:float).


BP5
~~~
Attention: BP5 in its original ACMG definition is not implemented.
Only the adaptation made in the gene-specific recommendations for BRCA1 and BRCA2 is implemented.

- **BP5_enigma**:
  Can be applied in case multifactorial likelihood analysis data is available for variant.
  Threshold for likelihood can be set under likelihood_thresholds for benignity and pathogenicity and differnet evidence strengths.

BP7
~~~
Check if deep intronic or synonymous variant does not affect splicing through a splicing prediction tool.
All implementation of BP7 check splice assay for benignity evidence.
If splice assay data is available results from the splicing assay are returned and prediction is not being assessed.

- **BP7**:
  Applies only to synonymous variants.

- **BP7_deep_intronic_ATM**:
  Applies to synonymous variants and deep intronic variants located at position >7 or <-40.

- **BP7_deep_intronic_enigma**:
  Applies to synonymous variants and deep intronic variants located at position >=7 or <=-21.
- **BP7_deep_intronic_enigma_check_disease_region**:
  Applies to synonymous variants and deep intronic variants located at position >=7 or <=-21.
  Checks that variant is located outside of coldspot region.
  Bed file with coldspot regions can be defined under annotation_files/critical_regions/coldspot_region.

- **BP7_deep_intronic_PALB2**:
  Applies to synonymous variants and deep intronic variants located at position >7 or <-21.
