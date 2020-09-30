#!/bin/bash
# This script concatenates two VCF files. It is important that they have the same samples and info structure

VCF1=$1; VCF2=$2; OUTPUT=$3

{ cat ${VCF1} | awk '!/^#/{exit}1'; cat ${VCF1} ${VCF2} | grep -v "^#" | sort -nk2,2; } > ${OUTPUT}
