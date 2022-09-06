## Tutorial for performing GWAS, according to requirements from CHOP

Please be aware that the following code has not been testing. If you have any questions or concerns, please do not hesitate to contact Matthew Patrick.

First, we set up the path to the phenotype and genotype data:
```
phenoDir=/nfs/turbo/umms-alextsoi/data/umich
genoDir=/nfs/turbo/precision-health/Genetic_Data/Freeze4/TOPMed/
```

### Extracting phenotype information

We need to get the patients who have particular ICD-9/10 codes. The parameter to this function is the name of a file with two columns (the ICD codes to use, and whether they are ICD9 or ICD10):
```
diag(){
  icdFile=$1
  zcat $phenoDir/diag.txt.gz | sed 1d \
  | awk -F'\t' 'NR==FNR{code[$1"_"$2]=1;next} \
    {pat[$1]} ($5"_"$6 in code){pat[$1]=1} \
    END{for(p in pat){print p"\t"(pat[p]==1)}' $icdFile -
}
```

We also need to get the patient demographics (here age and sex are used). The parameter to this function is self-identified race (for example "BLACK OR AFRICAN AMERICAN" or "WHITE OR CAUCASIAN"):
```
demo(){
  race=$1
  zcat demographics.txt.gz | awk -F'\t' -vrace=$race \
    '($4==race){print $1"\t"2021-$2"\t"$3}'
```

Then, we combine these two together, to get the phenotype/covariates for patients that have both:
```
pheno(){
  icdFile=$1; race=$2
  awk 'NR==FNR{diag[$1]=$2;next} ($1 in diag){print $0"\t"diag[$1]}' \
    <(diag $icdFile) <(demo $race)
}
```

The identifiers are encrypted, and we need to decrypt them using the crosswalk. However, please do not store any decrypted data, and delete it straight away:
```
7za e 'Crosswalk Data Direct DEID to Custom data DEID.zip' -p[password] -so \
| tr -d '"\r' | tr ',' '\t' | sed 1d | \
awk 'NR==FNR{id[$2]=$1;next} ($1 in id){print id[$1]"\t"$0}' \
  - <(pheno) | cut -f1,3- > .pheno
