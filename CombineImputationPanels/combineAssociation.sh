#!/bin/bash
MARKER_PATH=..
ASSOC_PATH=..
COHORTS="Cohort1 Cohort2"

panelMarkers(){
CHR=$1; PANEL=$2; COHORT=$3
awk -vpanel=${PANEL} 'NR==FNR&&$4==panel{markers[$1]=1;next} FNR>1&&markers[$3]' ${MARKER_PATH}/${COHORT}/data/markers${CHR}.dat ${ASSOC_PATH}/${COHORT}/results/plink2.${COHORT}_psv_v_ctl_${PANEL}_${CHR}.dat
}

for COHORT in $COHORTS; do
for CHR in $(seq 22); do
{ head -1 ${ASSOC_PATH}/${COHORT}/results/plink2.${COHORT}_psv_v_ctl_G1K_${CHR}.dat
{ panelMarkers ${CHR} G1K ${COHORT}
panelMarkers ${CHR} HRC ${COHORT}; } | sort -nk1,1 -nk2,2 -k4,4 -k5,5
} > ${ASSOC_PATH}/${COHORT}/results/plink2.${COHORT}_psv_v_ctl_${CHR}.dat
done
done
