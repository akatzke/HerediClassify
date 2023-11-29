# Automated variant classification

## Installation
Installation process has been tested for python 3.9 and 3.10.
Python dev version is needed.

1. Install general dependencies:
```sh
sudo apt install libpq-dev
```

2. Install python packages
```sh
pip install -r requirements.txt
```

3. Install pyensembl database
```sh
bash install_dependencies/pyensebmle_install.sh
```

4. Install non-python dependencies
The versions indicated have been tested with the tool
```sh
bash install_dependencies/install_bedtools.sh -p PATH -v 2.29.1
```
```sh
bash install_dependencies/install_htslib.sh -p PATH -v 1.18
```
```sh
bash install_dependencies/install_samtools.sh -p PATH -v 1.11
```

5. Download and format databases
The script will create a database folder with different subfolders per database.
The scripts expects the python dependencies installed above to be available.
```sh
bash install_dependencies/download_data.sh -p PATH
```

## Configuration
The file paths in the configuration file (config.yaml) need to be changed. Changing the root directory under annotation_files should suffice.

## Testing
Tests are implemented using pytest. To test general functionality execute pytest.

## Execution
```sh
python variant_classification/classify.py -c config.yaml -p json_string
```

## Execution via FastAPI
1. Start the FastAPI server
```sh
python webservice.py
```
This will start a uvicorn server running on http://0.0.0.0:8080
2. Execute classify on server
Send a curl request using the following command to the server
```sh
curl -X 'POST' \
  'http://0.0.0.0:8080/classify_variant' \
  -H 'accept: application/json' \
  -H 'Content-Type: application/json' \
  -d '{
  "config_path": "/home/katzkean/variant_classification/config.yaml",
  "variant_json": "{\"chr\": \"17\", \"pos\": 43057110, \"gene\": \"BRCA1\", \"ref\": \"A\", \"alt\": \"C\", \"variant_type\": [\"missense_variant\"], \"variant_effect\": [{\"transcript\": \"ENST00000357654\", \"hgvs_c\": \"c.5219T>G\", \"hgvs_p\": \"p.Val1740Gly\", \"variant_type\": [\"missense_variant\"], \"exon\": 19}, {\"transcript\": \"ENST00000471181\", \"hgvs_c\": \"c.5282T>G\", \"hgvs_p\": \"p.Val1761Gly\", \"variant_type\": [\"missense_variant\"], \"exon\": 20}], \"splicing_prediction_tools\": {\"SpliceAI\": 0.5}, \"pathogenicity_prediction_tools\": {\"REVEL\": 0.5, \"BayesDel\": 0.5}, \"gnomAD\": {\"AF\": 0.007, \"AC\": 12, \"popmax\": \"EAS\", \"popmax_AF\": 0.009, \"popmax_AC\": 5}, \"FLOSSIES\": {\"AFR\": 9, \"EUR\": 130}, \"mRNA_analysis\": {\"performed\": true, \"pathogenic\": true, \"benign\": true}, \"functional_data\": {\"performed\": true, \"pathogenic\": true, \"benign\": true}, \"prior\": 0.25, \"co-occurrence\": 0.56, \"segregation\": 0.56, \"multifactorial_log-likelihood\": 0.56, \"VUS_task_force_domain\": true, \"cancer_hotspot\": true, \"cold_spot\": true}"
}'
```
This will create the following output, with metadata included at the end
```sh
{
  "result": "{\"PVS1\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"very_strong\", \"comment\": \"PVS1 does not apply to this variant, as PVS1 does not apply to variant types [<VARTYPE.MISSENSE_VARIANT: 'missense_variant'>].\"}, \"PS1_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"strong\", \"comment\": \"No matches found for variant.\"}, \"PS1_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"strong\", \"comment\": \"No matches found for variant.\"}, \"PM1\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": true, \"strength\": \"moderate\", \"comment\": \"Variant in mutational hotspot.\"}, \"PM2\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"Variant occures with 0.007 in gnomAD.\"}, \"PM4\": {\"rule_type\": \"general\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"PM4 does not apply to this variant, as PVS1 does not apply to variant types [<VARTYPE.MISSENSE_VARIANT: 'missense_variant'>].\"}, \"PM5_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"No matches found for variant.\"}, \"PM5_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"moderate\", \"comment\": \"No matches found for variant.\"}, \"PP3_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"pathogenic\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Variant is not predicted to be pathogenic.\"}, \"PP3_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"pathogenic\", \"status\": true, \"strength\": \"supporting\", \"comment\": \"Variant is predicted to be pathogenic.\"}, \"BA1\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"stand_alone\", \"comment\": \"Variant occurs with 0.007 in gnomAD.\"}, \"BS1\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"strong\", \"comment\": \"Variant occurs with 0.007 in gnomAD.\"}, \"BS2\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": true, \"strength\": \"strong\", \"comment\": \"Something\"}, \"BP3\": {\"rule_type\": \"general\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"BP3 does not apply to this variant, as BP3 does not apply to variant types [<VARTYPE.MISSENSE_VARIANT: 'missense_variant'>].\"}, \"BP4_protein\": {\"rule_type\": \"protein\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Varinat is not predicted to be pathogenic.\"}, \"BP4_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Varinat is not predicted to be benign.\"}, \"BP7_splicing\": {\"rule_type\": \"splicing\", \"evidence_type\": \"benign\", \"status\": false, \"strength\": \"supporting\", \"comment\": \"Varinat is not predicted to be pathogenic.\"}}",
  "config": "config.yaml",
  "date": "2023-11-29",
  "version": "0.1.0"
}
```
