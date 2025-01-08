Execution
^^^^^^^^^^^^

HerediClassify can either be run in the command line, alternatively a FastAPI is implemented

Execution in the command line
==============================

.. code:: text
    python variant_classification/classify.py -c config.yaml -p json_string


Execution via FastAPI
======================

**1. Start the FastAPI server**

This will start a uvicorn server running on http://0.0.0.0:8080
.. code:: text
    python webservice.py

**2. Execute classify on server**

Send a curl request using the following command to the server
.. code:: text
    curl -X 'POST' \
    'http://0.0.0.0:8080/classify_variant' \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{
    "config_path": "/home/katzkean/variant_classification/config.yaml",
    "variant_json": "{\"chr\": \"17\", \"pos\": 43057110, \"gene\": \"BRCA1\", \"ref\": \"A\", \"alt\": \"C\", \"variant_type\": [\"missense_variant\"], \"variant_effect\": [{\"transcript\": \"ENST00000357654\", \"hgvs_c\": \"c.5219T>G\", \"hgvs_p\": \"p.Val1740Gly\", \"variant_type\": [\"missense_variant\"], \"exon\": 19}, {\"transcript\": \"ENST00000471181\", \"hgvs_c\": \"c.5282T>G\", \"hgvs_p\": \"p.Val1761Gly\", \"variant_type\": [\"missense_variant\"], \"exon\": 20}], \"splicing_prediction_tools\": {\"SpliceAI\": 0.5}, \"pathogenicity_prediction_tools\": {\"REVEL\": 0.5, \"BayesDel\": 0.5}, \"gnomAD\": {\"AF\": 0.007, \"AC\": 12, \"popmax\": \"EAS\", \"popmax_AF\": 0.009, \"popmax_AC\": 5}, \"FLOSSIES\": {\"AFR\": 9, \"EUR\": 130}, \"mRNA_analysis\": {\"performed\": true, \"pathogenic\": true, \"benign\": true}, \"functional_data\": {\"performed\": true, \"pathogenic\": true, \"benign\": true}, \"prior\": 0.25, \"co-occurrence\": 0.56, \"segregation\": 0.56, \"multifactorial_log-likelihood\": 0.56, \"VUS_task_force_domain\": true, \"cancer_hotspot\": true, \"cold_spot\": true}"
    }'
This will create the following output, with metadata included at the end

.. code:: text
    {"result":"{\"PVS1\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"very_strong\", \"comment\": \"PVS1 does not apply to this variant, as PVS1 does not apply to variant types missense_variant.\"}, \"PS1_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"strong\", \"comment\": \"No ClinVar entries found that show the same amino acid change as pathogneic.\"}, \"PS1_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"strong\", \"comment\": \"No ClinVar entries found that show splice variants at the same nucleotide position as pathogenic..\"}, \"PM1\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": true, \"strength\": \"moderate\", \"comment\": \"Variant in mutational hotspot.\"}, \"PM2\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"Variant occures with 0.009 in GnomAD subpopulation EAS.\"}, \"PM4\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"PM4 does not apply to this variant, as PVS1 does not apply to variant types missense_variant.\"}, \"PM5_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"No ClinVar entries found that show an amino acid change in the same position as pathogenic.\"}, \"PM5_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"No ClinVar entries found that show variant in the same splice site as pathogenic.\"}, \"PP3_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Variant is not predicted to be pathogenic by REVEL.\"}, \"PP3_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"pathogenic\", \"status\": true, \"strength\": \"supporting\", \"comment\": \"Variant is predicted to have a splice effect by SpliceAI.\"}, \"BA1\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"stand_alone\", \"comment\": \"Variant occures with 0.009 in GnomAD subpopulation EAS.\"}, \"BS1\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"strong\", \"comment\": \"Variant occures with 0.009 in GnomAD subpopulation EAS.\"}, \"BS2\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": true, \"strength\": \"strong\", \"comment\": \"The variant occures 130 in FLOSSIES.\"}, \"BP3\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"BP3 does not apply to this variant, as BP3 does not apply to variant types missense_variant.\"}, \"BP4_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Variant is not predicted to be benign REVEL.\"}, \"BP4_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Variant is not predicted to have no splicing effect by SpliceAI.\"}, \"BP7_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Variant is not predicted to have no splicing effect by SpliceAI.\"}}",
    "config_file":"config.yaml",
    "config_name":"acmg brca2",
    "date":"2023-12-04",
    "version":"0.1.0"}
