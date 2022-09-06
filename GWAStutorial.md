## Tutorial for performing GWAS, according to requirements from CHOP

Please be aware that the following code has not been testing. If you have any questions or concerns, please do not hesitate to contact Matthew Patrick.

### Extracting phenotype information

```
diag(){
icdFile=$1
zcat diag.txt.gz | sed 1d \
| awk -F'\t' 'NR==FNR{code[$1"_"$2]=1;next} \
{pat[$1]} ($5"_"$6 in code){pat[$1]=1} \
END{for(p in pat){print p"\t"(pat[p]==1)}' $icdFile -
}
```
