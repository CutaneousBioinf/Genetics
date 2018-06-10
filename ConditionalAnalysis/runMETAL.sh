#!/bin/bash

STEP=$1

/usr/cluster/bin/metal << EOT

SCHEME STDERR
STDERRLABEL STDERR

MARKER MARKER
EFFECT BETA
STDERR  SE

GENOMICCONTROL 1.025
PROCESS ../results/${STEP}/PsAGWAS.psa_v_ctl.dat

GENOMICCONTROL 1.015
PROCESS ../results/${STEP}/CASP.psa_v_ctl.dat

GENOMICCONTROL 1.045
PROCESS ../results/${STEP}/KielGWAS.psa_v_ctl.dat

GENOMICCONTROL 1.015
PROCESS ../results/${STEP}/Genizon.psa_v_ctl.dat

GENOMICCONTROL 1.02
PROCESS ../results/${STEP}/Exomechip.psa_v_ctl.dat

OUTFILE ../results/${STEP}/meta.psa_v_ctl. tbl
ANALYZE HETEROGENEITY
QUIT
EOT

/usr/cluster/bin/metal << EOT

SCHEME STDERR
STDERRLABEL STDERR

MARKER MARKER
EFFECT BETA
STDERR  SE

GENOMICCONTROL 1.015
PROCESS ../results/${STEP}/CASP.psc_v_ctl.dat

GENOMICCONTROL 1.045
PROCESS ../results/${STEP}/KielGWAS.psc_v_ctl.dat

GENOMICCONTROL 1.015
PROCESS ../results/${STEP}/Genizon.psc_v_ctl.dat

GENOMICCONTROL 1.02
PROCESS ../results/${STEP}/Exomechip.psc_v_ctl.dat

OUTFILE ../results/${STEP}/meta.psc_v_ctl. tbl
ANALYZE HETEROGENEITY
QUIT
EOT

Rscript <(join -j1 <(cat ../results/${STEP}/meta.psa_v_ctl.1tbl | sed 1d | cut -f1,4,5 | sort -k1,1) <(cat ../results/${STEP}/meta.psc_v_ctl.1tbl | sed 1d | cut -f1,4,5 | sort -k1,1) | awk '{print "cat(\"" $1 "\",pchisq(" ($2-$4)^2/($3^2+$5^2) ",df=1,lower.tail=F),\"\\n\")" }') | grep -v "NA" | sort -gk2,2 > ../results/${STEP}/PsAPsC.indirect
