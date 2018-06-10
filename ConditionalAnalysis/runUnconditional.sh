#!/bin/bash
function runAssociation(){
COHORT=$1; TEST=$2; STEP=$3
plink2 --pfile ../data/${COHORT}_${TEST} --pheno ../data/${COHORT}_${TEST}.phe --pheno-name ${TEST} --covar ../data/${COHORT}_${TEST}.covar --glm hide-covar cols=+beta --threads 10 --out ../results/${STEP}/plink2.${COHORT}
awk 'BEGIN{print "MARKER\tBeta\tSE"} NR>1{print $3"\t"$8"\t"$9}' ../results/${STEP}/plink2.${COHORT}.${TEST}.glm.logistic | sort -nk1,1 > ../results/${STEP}/${COHORT}.${TEST}.dat
}

STEP=0
mkdir ../results/${STEP}
for COHORT in CASP Exomechip Genizon KielGWAS PsAGWAS; do
runAssociation ${COHORT} psa_v_ctl ${STEP}; done
for COHORT in CASP Exomechip Genizon KielGWAS; do
runAssociation ${COHORT} psc_v_ctl ${STEP}; done
./runMETAL.sh ${STEP}
BESTMARKER=$(cat ../results/${STEP}/PsAPsC.indirect | awk 'BEGIN{getline;print $1}')
echo ${BESTMARKER} > ../data/condition.dat
