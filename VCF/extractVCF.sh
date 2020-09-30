#!/bin/bash
# This script extracts specific markers from a VCF file to create a new VCF subset. The marker format must be CHROM:POS_REF/ALT_ID

VCF=$1; MARKERS=$2; OUTPUT=$3

extractHeader(){ cat $1 | awk '!/^#/{exit}1'; }
stripHeader(){ cat $1 | grep -v "^#"; }
extractMarkers(){ awk 'FNR==NR{markers[$0]=1;next} markers[$1":"$2"_"$4"/"$5"_"$3]{print $0}' $1 -; }

{ extractHeader ${VCF}; stripHeader ${VCF} | extractMarkers ${MARKERS}; } > ${OUTPUT}
