#!/bin/bash
IMP_PATH=..
OUT_PATH=..
COHORTS="Cohort1 Cohort2"

imputationQuality(){
CHR=$1; PANEL=$2; COHORT=$3
awk 'NR>1&&$7>=.7{print $1"_"$2"/"$3"_"$1"\t"$7}' ${IMP_PATH}/${COHORT}/${PANEL}/${COHORT}_chr${CHR}.shapeIT2.minimac3_${PANEL}.info | sort -k1,1
}

for COHORT in $COHORTS; do
for CHR in $(seq 22); do
join -j1 -a1 -a2 -e0 -o0,1.2,2.2 <(imputationQuality ${CHR} G1K ${COHORT}) <(imputationQuality ${CHR} HRC ${COHORT}) | awk '$2>$3{print $0"\tG1K";next} {print $0"\tHRC"}' > ${OUT_PATH}/${COHORT}/data/markers${CHR}.dat
done
done
