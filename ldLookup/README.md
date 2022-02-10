
# ldLookup

This software facilitates efficient lookup and analysis of linkage disequilibrium (LD) between genetic variants.

## Installation
To install, clone this repo with `git clone  https://github.com/CutaneousBioinf/Genetics.git`, change directory to `ldLookup`, and run `make`. This will create the `ldLookup` executable, which you can copy as needed.

```
>>> git clone https://github.com/CutaneousBioinf/Genetics.git
>>> cd Genetics/ldLookup
>>> make
...
>>> ls
... ldLookup ...
```

## Usage
This section shows how to use `ldLookup`. `ldLookup` provides context-aware help messages via the `-h` or `--help` flags.
```
>>> ./ldLookup --help
ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants
Usage: ./ldLookup [OPTIONS] SUBCOMMAND

Options:
  -h,--help                   Print this help message and exit

Subcommands:
  create                      Create a new dataset from LD data
  get_ld                      Get markers in LD with a particular index marker
  similar_by_value            Retrieve markers with MAF and number of LD surrogates near target values
  similar_by_snp              Retrieve markers with MAF and number of LD surrogates similar to a key marker
  distribute                  Generate distributions for hypothesis testing

>>> ./ldLookup create -h
Create a new dataset from LD data
Usage: ./ldLookup create [OPTIONS] dir path_to_data

Positionals:
  dir TEXT REQUIRED           Specifies a directory that will contain the created dataset
  path_to_data TEXT REQUIRED  File containing LD data

Options:
  -h,--help                   Print this help message and exit
  --dir TEXT REQUIRED         Specifies a directory that will contain the created dataset
  -p,--path-to-data TEXT REQUIRED
                              File containing LD data
  -t,--threshold FLOAT=0      LD score (as measured by R^2) below which a marker will be ignored
  -I,--index-snp-id-column UINT=3
                              Index to column of data containing index SNP IDs
  -M,--maf-col UINT=4         Index to column of data containing MAFs for index SNPs
  -L,--ld-snp-id-column UINT=7
                              Index to column of data containing IDs for SNPs in LD with the index SNP
  -R,--r2-index UINT=9        Index to column of data containing r-squared values
  -d,--delimiter CHAR=        Character used to separate columns of data
  -k,--key-size-limit UINT=200
                              Maximum index SNP ID length in bytes
  -l,--ld-bin-count UINT=15   Approximate number of groupings for index SNPs by number of LD surrogates
  -m,--maf-bin-count UINT=15  Approximate number of groupings for index SNPs by MAF

>>> ./ldLookup get_ld -h
Get markers in LD with a particular index marker
Usage: ./ldLookup get_ld [OPTIONS] dir [path_to_markers] [markers...]

Positionals:
  dir TEXT REQUIRED           Specifies the location of the dataset to operate on (see 'ldLookup create')
  path_to_markers TEXT        File containing newline-separated index SNP IDs
  markers TEXT=[] ...         Additional space-separated SNP IDs

Options:
  -h,--help                   Print this help message and exit
  --dir TEXT REQUIRED         Specifies the location of the dataset to operate on (see 'ldLookup create')
  -p,--path-to-markers TEXT   File containing newline-separated index SNP IDs
  -m,--markers TEXT=[] ...    Additional space-separated SNP IDs

>>>./ldLookup similar_by_value -h
Retrieve markers with MAF and number of LD surrogates near target values
Usage: ./ldLookup similar_by_value [OPTIONS] dir target_maf target_surrogates

Positionals:
  dir TEXT REQUIRED           Specifies the location of the dataset to operate on (see 'ldLookup create')
  target_maf FLOAT=0 REQUIRED Target MAF value
  target_surrogates UINT=0 REQUIRED
                              Target number of LD surrogates

Options:
  -h,--help                   Print this help message and exit
  --dir TEXT REQUIRED         Specifies the location of the dataset to operate on (see 'ldLookup create')
  -m,--target-maf FLOAT=0 REQUIRED
                              Target MAF value
  -s,--target-surrogates UINT=0 REQUIRED
                              Target number of LD surrogates

>>> ./ldLookup similar_by_snp -h
Retrieve markers with MAF and number of LD surrogates similar to a key marker
Usage: ./ldLookup similar_by_snp [OPTIONS] dir [path_to_markers] [markers...]

Positionals:
  dir TEXT REQUIRED           Specifies the location of the dataset to operate on (see 'ldLookup create')
  path_to_markers TEXT        File containing newline-separated index SNP IDs
  markers TEXT=[] ...         Additional space-separated SNP IDs

Options:
  -h,--help                   Print this help message and exit
  --dir TEXT REQUIRED         Specifies the location of the dataset to operate on (see 'ldLookup create')
  -p,--path-to-markers TEXT   File containing newline-separated index SNP IDs
  -m,--markers TEXT=[] ...    Additional space-separated SNP IDs

>>> ./ldLookup distribute -h
Generate distributions for hypothesis testing
Usage: ./ldLookup distribute [OPTIONS] dir [path_to_markers] [markers...]

Positionals:
  dir TEXT REQUIRED           Specifies the location of the dataset to operate on (see 'ldLookup create')
  path_to_markers TEXT        File containing newline-separated index SNP IDs
  markers TEXT=[] ...         Additional space-separated SNP IDs

Options:
  -h,--help                   Print this help message and exit
  --dir TEXT REQUIRED         Specifies the location of the dataset to operate on (see 'ldLookup create')
  -p,--path-to-markers TEXT   File containing newline-separated index SNP IDs
  -m,--markers TEXT=[] ...    Additional space-separated SNP IDs
  -d,--distributions UINT=1   Number of distributions to generate
```

### A Brief Example
```
>>> head tst/test.ld
CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2 
1 11008 1:11008:C:G 0.0884692 1 11012 1:11012:C:G 0.0884692 1 
1 46285 1:46285:ATAT:A 0.000994036 1 66461 1:66461:T:A 0.000994036 1 
1 46285 1:46285:ATAT:A 0.000994036 1 81590 rs202072409:81590:AC:A 0.000994036 1 
1 13116 rs201725126:13116:T:G 0.186879 1 13118 rs200579949:13118:A:G 0.186879 1 
1 14599 1:14599:T:A 0.161034 1 14604 1:14604:A:G 0.161034 1 
1 15274 rs201931625:15274:A:G 0.292247 1 15274 rs201931625:15274:A:T 0.292247 1 
1 15644 1:15644:G:A 0.00795229 1 15774 rs374029747:15774:G:A 0.00795229 1 
1 17571 1:17571:C:T 0.000994036 1 49343 1:49343:T:C 0.000994036 1 
1 49554 1:49554:A:G 0.0636183 1 76838 1:76838:T:G 0.0616302 0.889466 
>>> ./ldLookup create my_dataset tst/test.ld
Did not parse malformed line CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2 
Did not parse malformed line CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2 
Did not parse malformed line CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2
>>> ./ldLookup get_ld my_dataset -m 1:11008:C:G 1:46285:ATAT:A
1:11008:C:G     1:11012:C:G
1:46285:ATAT:A  1:66461:T:A
1:46285:ATAT:A  rs202072409:81590:AC:A
>>> ./ldLookup get_ld my_dataset -m bad_key
ERROR: VectorDiskHash Error in 'ld_table': Key Error - Key 'bad_key' does not exist
>>>./ldLookup similar_by_value my_dataset --target-maf 0.08 --target-surrogates 1
0.080000 1      1:11008:C:G
0.080000 1      1:60351:A:G
0.080000 1      rs184538873:77874:G:A
0.080000 1      rs375955515:89599:A:T
0.080000 1      rs147538909:115746:C:T
0.080000 1      rs7349153:565490:T:C
0.080000 1      rs201234755:701779:GAATA:G
0.080000 1      rs116587930:727841:G:A
0.080000 1      rs12131618:732809:T:C
0.080000 1      rs376110660:804188:G:C
0.080000 1      rs72631880:805556:T:A
>>> ./ldLookup similar_by_snp my_dataset -m 1:11008:C:G
1:11008:C:G     1:11008:C:G
1:11008:C:G     1:60351:A:G
1:11008:C:G     rs184538873:77874:G:A
1:11008:C:G     rs375955515:89599:A:T
1:11008:C:G     rs147538909:115746:C:T
1:11008:C:G     rs7349153:565490:T:C
1:11008:C:G     rs201234755:701779:GAATA:G
1:11008:C:G     rs116587930:727841:G:A
1:11008:C:G     rs12131618:732809:T:C
1:11008:C:G     rs376110660:804188:G:C
1:11008:C:G     rs72631880:805556:T:A
>>> ./ldLookup distribute my_dataset -m 1:11008:C:G -d 5
Distribution #1 rs116587930:727841:G:A
Distribution #2 rs12131618:732809:T:C
Distribution #3 rs201234755:701779:GAATA:G
Distribution #4 rs376110660:804188:G:C
Distribution #5 1:11008:C:G
```

### Create
`ldLookup` has five subcommands. The first one to use is `create`, which uses files of linkage disequilibrium data to build a dataset `ldLookup` can use for sampling. This will not modify the original data (although it's always advisable to make a back-up). For example, this command builds a dataset in the directory `example_dataset` containing LD data from `data/example.ld`.

```
>>> ./ldLookup create example_dataset data/example.ld 
>>> ll
...
example_dataset
...
```
`create` supports different file formats via additional options. Run `./ldLookup create --help` for a complete list.
- The flags `-I`, `-M`, `-L`, and `-R` support different column orderings. They are case-sensitive.
- The flag `-d` supports different delimiters between columns. Note that the software treats consecutive delimiters as one delimiter, which works as expected for columns separated by different amounts of whitespace. It does not work as expected for `.csv` files since `A,B,C` and `A,B,,C` are considered identical. 
- The flag `-k` changes the maximum size of data in the `INDEX_SNP` column. If you get an error like `Key Too Large`, increase `-k`.
- The flag `-t` sets a minimum value for r-squared. If some row of the input data has an r-squared value below this minimum, that row will not be part of the output dataset. This defaults to 0, so it should probably be set.
- The flags `-l` and `-m` control sampling granularity. Higher values are more granular.

Two datasets cannot be assigned the same name. So, running ``./ldLookup create example_dataset`` a second time will produce an error. To free a name, manually delete the dataset's files (e.g. ``rm example_dataset_*``).

### distribute
This subcommand controls `ldLookup`'s main function: randomly generating empirical distributions while controlling for MAF and number of LD surrogates.

This subcommand is not yet documented (beyond the help messages)

### get_ld / similar_by_snp / similar_by_value
These subcommands are not yet documented (beyond the help messages).
