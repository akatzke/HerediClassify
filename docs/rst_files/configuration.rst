Configuration
^^^^^^^^^^^^^^^^^^^^^^^^

A number of different thresholds and paths can be defined in the configuration file.
An example configuration ~config.yaml~ is provided in the repository.
Additionally, configuration files for the gene specific recommendations for the hereditary breast and ovarian cancer risk genes ATM (v.1.3.0), BRCA1 (v.1.1.0), BRCA2 (v.1.1.0), CDH1 (v.3.1.0), PALB2 (v.1.1.0), PTEN (v.3.1.0) and TP53 (v.1.4.0) are provided (see ~gene_specific~ folder).

Configuration name
====================
At the top at each configuration file a name and a version can be defined.
Version and name of the configuration together with the final classification result.

Rules
=======
In the rules section of the configuration file, a list of rules from the implemented rules can be defined in a list.
The processing of the rules is case insensitive.

Disease relevant transcript
=============================
In this section of the configuration file a disease relevant transcript or list of disease relevant transcripts can be defined.
The ~name~ needs to be the Ensembl transcript ID without version extension.
Under ~nmd_threshold~ the last nucleotide position that may cause nonsense mediated decay can be documented.
~pos_last_known_patho_ptc~ allows for the definition of the last known pathogenic premature termination codon.
These transcript specific information are often needed for the application of PVS1.

Thresholds prediction tools
========================================
In this section thresholds for pathogenicitiy prediction tools can be defined.
HerediClassify is only designed to use on prediction tool for pathogenicity prediction and one prediction tool for splicing prediction.
For pathogenicity prediction the thresholds can be defined under ~pathogenicity_prediction~ and thresholds for splicing prediction tools can be defined under ~splicing prediction~.
Both thresholds have the same structure.
Under ~name~ the name of the prediction tool can be set.
The prediction tool should be the same as the one specific in the variant import json.
Under ~benign~ and ~pathogenic~ the thresholds for pathogenic and benign evidence can be set.
Using ~direction~ the direction of the comparison between the threshold and the prediction value is defined.
The following directions are possible:
- ~less~
- ~less_than_or_equal~
- ~greater~
- ~greater_than_or_equal~
To define the thresholds use evidence strength: threshold e.g. ~supporting:0.1~.
When using the rules ~PP3_protein_mult_strength~ or ~BP4_protein_mult_strength~ thresholds for higher evidence strengths can be defined using ~moderate~, ~strong~ and ~very_strong~.
For example see example configuration files under ~gene_sepcific/pejaver_mult_strength~.
When another implementation of PP3 or BP4 is used that does not apply an evidence strength higher than supporting, thresholds for higher evidence strengths will be ignored during the assessment of the criterion.

Likelihood thresholds
========================================
Under likelihood thresholds the threshold for all likelihoods can be set.
These can be accessed by various criteria including BS4 and PP1.

Allele frequency thresholds
========================================
Here allele frequency thresholds can be defined as decimal figures.
CAVE no conversion into percentage is implemented.
- **threshold_ba1**
  Threshold for BA1.
- **threshold_bs1**
  Threshold for BS1.
- **threshold_bs1_supporting**
  Threshold for BS1 with supporting evidence, required for BS1_supporting
- **threshold_bs2**
- **threshold_bs2_supporting**
- **threshold_pm2**
  Threshold for PM2.
