## Tutorial for performing GWAS, according to requirements from CHOP

Please be aware that the following code has not been testing. If you have any questions or concerns, please do not hesitate to contact Matthew Patrick.

First, we set up the path to the phenotype and genotype data:
```
phenoDir=/nfs/turbo/umms-alextsoi/data/umich
genoDir=/nfs/turbo/precision-health/Genetic_Data/Freeze4
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
  zcat $phenoDir/demographics.txt.gz | awk -F'\t' -vrace=$race \
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
  - <(pheno) | cut -f1,3- | awk 'BEGIN{print "FID\tIID\tage\tsex\tdis"} {print "0\t"$0}' > dis.pheno
```

We extract the genotyped data for analysis:
```
plink2 --bfile $genoDir/directly-typed/MGI.Freeze4.DIRECTLY.TYPED.hg38.bim --const-fid --keep dis.pheno --make-bed --out typed --maf 0.01 --geno 0.02 --mind 0.02
```

PCA is performed to check whether the genetic ancestry for any patients differs substantially from the self-identified ancestry:
```
plink-1.9 --bfile typed --indep-pairwise 100 5 0.2
plink-1.9 --bfile typed --extract plink.prune.in --make-bed --out typed_LD
plink-1.9 --bfile typed_LD --pca --out typed.pc
awk -vOFS='\t' 'NR==FNR{pc[$2]=$3"\t"$4"\t"$5;next} FNR==1{$NF="PC1\tPC2\tPC3";print;next} {print "0\t"$1"\t"$2"\t"$3"\t"pc[$1]}' typed.pc dis.pheno > dis.covar
```

Next, we extract the imputed genetic data for these patients:
```
for chr in {1..22}; do
  plink2 --vcf $genoDir/TOPMed/TOPMED.FILTERED.chr$chr.vcf.gz dosage=DS --const-fid --keep dis.pheno --make-pgen --out chr$chr --maf 0.01
done
```

Using plink, we perform association analysis:
```
for chr in {1..22}; do
  plink2 --pfile chr$chr --pheno dis.phe --pheno-name dis --covar dis.covar --glm hide-covar cols=+beta --threads 10 --memory 1024 --out dis
done
```

Using REGENIE, we perform association analysis:
```
regenie --step 1 --bed typed --covarFile dis.covar --phenoFile dis.phe --bsize 100 --bt --lowmem --lowmem-prefix tmp_rg --out dis
regenie --step 2 --pgen dis.pgen --covarFile dis.covar --phenoFile dis.phe --bsize 200 --bt --firth --approx --pThresh 0.01 --pred dis_pred.list --out dis
```
