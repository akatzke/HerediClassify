#!/bin/bash
# Install pyensemble refernce genome
set -e
set -o pipefail


helpFunction()
{
   echo ""
   echo "Usage: $0 -v version"
   echo "This script installs Ensembl release data"
   echo -e "\t-v The ensembl release version. Eg. 110"
   exit 1 # Exit script after printing help
}

while getopts "v:" opt
do
   case "$opt" in
      v ) version="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$version" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "Installing Ensembl release $version ..."

pyensembl install --release $version --species human
