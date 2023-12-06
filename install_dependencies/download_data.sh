#!/bin/bash
# Copied from GC-HBOC/HerediVar install scripts
set -e
set -o pipefail

helpFunction()
{
   echo ""
   echo "Usage: $0 -p path"
   echo "This script downloadd data"
   echo -e "\t-p The database directory"
   exit 1 # Exit script after printing help
}

while getopts "p:" opt
do
   case "$opt" in
      p ) basedir="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$basedir" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

dbs=$basedir/databases
mkdir -p $dbs

# Begin script in case all parameters are correct
echo "Downloading data to $dbs..."


cd $dbs
mkdir -p Clinvar
cd Clinvar
clinvar=$dbs/Clinvar

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
in_path_clinvar=$clinvar/clinvar.vcf.gz
python $basedir/variant_classification/install_dependencies/data_filter_clinvar.py -i $in_path_clinvar

cd $dbs
mkdir -p Uniprot
cd Uniprot
unip=$dbs/Uniprot


curl -L "http://api.genome.ucsc.edu/getData/track?genome=hg38;track=unipRepeat" > repeats_hg38_uniprot.json
in_path_uniprot=$unip/repeats_hg38_uniprot.json
python $basedir/variant_classification/install_dependencies/data_format_uniprot_rep.py -i $in_path_uniprot
