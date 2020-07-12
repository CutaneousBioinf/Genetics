#!/bin/bash
MARKER_PATH=..
VCF_PATH=..
COHORTS="Cohort1 Cohort2"

combineVCF(){ COHORT=$1; CHR=$2
{ tabix -H ${VCF_PATH}/${COHORT}/G1K/${COHORT}_chr22.shapeIT2.minimac3_G1K.dose.vcf.r7.gz
awk '{split($1,a,"_");split(a[2],b,"/");split(a[1],a,":");print a[1]"\t"a[2]"\t"b[1]"\t"b[2]"\t"$4}' ${MARKER_PATH}/${COHORT}/data/markers${CHR}.dat | \
while read CHR POS A1 A2 PANEL; do
   tabix ${VCF_PATH}/${COHORT}/${PANEL}/${COHORT}_chr${CHR}.shapeIT2.minimac3_${PANEL}.dose.vcf.r7.gz ${CHR}:${POS}-${POS} | awk -va1=${A1} -va2=${A2} '$4==a1&&$5==a2{print;exit}'
done
} | bgzip -c > ${VCF_PATH}/${COHORT}/merged/${COHORT}_chr${CHR}.shapeIT2.minimac3_merged.dose.vcf.r7.gz
tabix -pvcf ${VCF_PATH}/${COHORT}/merged/${COHORT}_chr${CHR}.shapeIT2.minimac3_merged.dose.vcf.r7.gz
}

for COHORT in $COHORTS; do
for CHR in $(seq 22); do
combineVCF ${COHORT} ${CHR}
done; done
