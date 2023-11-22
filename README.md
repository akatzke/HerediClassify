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
