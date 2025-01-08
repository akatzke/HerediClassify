Installation
^^^^^^^^^^^^

Installation process has been tested for python 3.9 and 3.10.
Python dev version is needed.

Install HerediClassify
========================

1. Install general dependencies:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sudo apt install libpg-dev

2. Install python packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    pip install -r requirements.txt

3. Install pyensembl database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   bash install_dependencies/install_pyensembl.sh

4. Install non-python dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The versions indicated have been tested with the tool

.. code:: bash

    bash install_dependencies/install_bedtools.sh -p PATH -v 2.29.1
    bash install_dependencies/install_htslib.sh -p PATH -v 1.18
    bash install_dependencies/install_samtools.sh -p PATH -v 1.11

Download additional data
========================

1. Download and format databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The script will create a database folder with different subfolders per database.
The scripts expects the python dependencies installed above to be available.

.. code:: bash

    bash install_dependencies/download_data.sh -p PATH

2. Annotate ClinVar with SpliceAI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A ClinVar file annotated with SpliceAI scores is needed for the application of PS1 and PM5.
In order to obtain this file manual steps and additional dependencies need to be resolved.

**1. Download SpliceAI scores from Illumina**
The SpliceAI scores can be found [here](https://basespace.illumina.com/analyses/194103939/files?projectId=66029966).
Please download the masked indel ans snv file for GRCh38:
- spliceai_scores.masked.indel.hg38.vcf.gz
- spliceai_scores.masked.indel.hg38.vcf.gz.tbi
- spliceai_scores.masked.snv.hg38.vcf.gz
- spliceai_scores.masked.snv.hg38.vcf.gz.tbi

**2. Install ngs-bits**
See the ngs-bits [github](https://github.com/imgag/ngs-bits) for install instructions.

**3. Perform annotation**
Use the merge_clinvar_spliceai.sh script to annotate the ClinVar download with the SpliceAI annotations. Make sure to change the paths at the top of the script to match paths on your system.

**4. Filtering**
Filter the ClinVar file annotated with SpliceAI by using the data_filter_clinvar.py script.

.. code:: bash

    python ../HerediClassify/install_dependencies/data_filter_clinvar.py -i path_to/clinvar_spliceai_all_sorted.vcf.gz

Testing
========
Tests are implemented using pytest. To test general functionality execute:

.. code:: bash

    pytest test/
