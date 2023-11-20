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
bash install_dependencies/install_htslib.sh -p PATH -v 
```
```sh
bash install_dependencies/install_samtools.sh -p PATH -v
```

5. Download and format databases
The script will create a database folder with different subfolders per database.
The scripts expects the python dependencies installed above to be available.
```sh
bash install_dependencies/download_data.sh -p PATH
```

## Configuration
The file paths in the configuration file (config.yaml) need to be changed for the new install.
Under annotation_files the root directory for the database folders created under 5. in the install instructions.

## Execution
```sh
python bin/classify.py -c config.yaml -p 
```
