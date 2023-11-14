#!/usr/bin/env sh

# Download json schema for
curl -L "http://api.genome.ucsc.edu/list/schema?genome=hg38;track=unipRepeat" > /home/katzkean/Downloads/repeats_schema.json

# Download file from USCS as json
curl -L "http://api.genome.ucsc.edu/getData/track?genome=hg38;track=unipRepeat" > /home/katzkean/Downloads/repeats_hg38_uniprot.json

# Here needs to be script that converts the Download to a bed file usable by me
