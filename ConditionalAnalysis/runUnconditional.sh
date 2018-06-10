#!/bin/bash
FOLD=$1
function runAssociation(){
COHORT=$1; TEST=$2; STEP=$3
plink2 --pfile ../data/${COHORT}_${TEST}${FOLD} --pheno ../data/${COHORT}_${TEST}_${FOLD}.phe --pheno-name ${TEST} --covar ../data/${COHORT}_${TEST}_${FOLD}.covar --glm hide-covar cols=+beta --threads 10 --out ../results/${STEP}/plink2.${COHORT}${FOLD}
awk 'BEGIN{print "MARKER\tBeta\tSE"} NR>1{print $3"\t"$8"\t"$9}' ../results/${STEP}/plink2.${COHORT}${FOLD}.${TEST}.glm.logistic | sort -nk1,1 > ../results/${STEP}/${COHORT}.${TEST}${FOLD}.dat
}

STEP=0
mkdir ../results/${STEP}
for COHORT in CASP Exomechip Genizon KielGWAS PsAGWAS; do
runAssociation ${COHORT} psa_v_ctl ${STEP}; done
for COHORT in CASP Exomechip Genizon KielGWAS; do
runAssociation ${COHORT} psc_v_ctl ${STEP}; done
./runMETAL.sh ${STEP} ${FOLD}
BESTMARKER=$(cat ../results/${STEP}/PsAPsC${FOLD}.indirect | awk 'BEGIN{getline;print $1}')
echo ${BESTMARKER} > ../data/condition${FOLD}.dat
