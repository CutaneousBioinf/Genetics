#!/bin/bash
# This script merges two VCF files (i.e. combines their samples). It is important that they have the same info structure

VCF1=$1; VCF2=$2; OUTPUT=$3

getHeader(){ cat $1 | awk '!/^#/{exit}1' | head -n -1; }
stripHeader(){ cat $1 | grep -vP '^#(?!CHROM)'; }

{ getHeader ${VCF1}; paste <(stripHeader ${VCF1}) <(stripHeader ${VCF2} | cut -f10-); } > ${OUTPUT}
