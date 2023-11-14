#!/bin/bash
# Fashioned after GC-HBOC/HerediVar install scripts
set -e
set -o pipefail


helpFunction()
{
    echo ""
    echo "Usage: $0 -p path -v version"
    echo "This script installs bedtools"
    echo -e "\t-p The path bedtools will be installed"
    echo -e "\t-v The bedtools version. E.g. 2.29.1"
   exit 1 # Exit script after printing help
}

while getopts "p:v:" opt
do
    case "$opt" in
        p ) basedir="$OPTARG" ;;
        v ) version="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "$basedir" ] || [ -z "$version" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "Installing bedtools $version to $basedir..."



cd $basedir
wget -q https://github.com/arq5x/bedtools2/releases/download/v$version/bedtools-$version.tar.gz
tar -zxvf bedtools-$version.tar.gz
cd bedtools2
make
