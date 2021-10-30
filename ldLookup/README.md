# ldLookup

This software facilitates efficient lookup and analysis of linkage disequilibrium (LD) between genetic variants.

## Installation
To install, clone this repo with `git clone  https://github.com/CutaneousBioinf/Genetics.git`, change directory to `ldLookup`, and run `make`. This will create the `ldlookupx` executable.

```
>>> git clone https://github.com/CutaneousBioinf/Genetics.git
>>> cd Genetics/ldLookup
>>> make
...
>>> ls
... ldLookupx ...
```

## Usage
```
>>> ./LdLookupx --help
ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants
Usage: ./ldLookupx [OPTIONS] table SUBCOMMAND

Positionals:
  table TEXT REQUIRED         Name of lookup table to search or create

Options:
  -h,--help                   Print this help message and exit
  -t,--t TEXT REQUIRED        Name of lookup table to search or create

Subcommands:
  create                      Create new lookup table
  retrieve                    Retrieve values existing from lookup table

>>> ./ldLookupx table create --help
Create new lookup table
Usage: ./ldLookupx create [OPTIONS] source_path

Positionals:
  source_path TEXT REQUIRED   File containing source data

Options:
  -h,--help                   Print this help message and exit
  -r,--r2-threshold FLOAT     Minimum r-squared value to include a key-value pair in the table
  -K,--key-index UINT         Zero-based index to column of data containing lookup table keys
  -V,--value-index UINT       Zero-based index to column of data containing lookup table values
  -R,--r2-index UINT          Zero-based index to column of data containing r-squared values
  -d,--delimiter CHAR         Character used to separate columns of data

>>> ./ldLookupx table retrieve --help
Retrieve values existing from lookup table
Usage: ./ldLookupx retrieve [OPTIONS] [file] [keys...]

Positionals:
  file TEXT                   File containing alleles to look up (one per line)
  keys TEXT ...               Alleles to look up

Options:
  -h,--help                   Print this help message and exit
  -f,--file TEXT              File containing alleles to look up (one per line)
  -k,--keys TEXT ...          Alleles to look up
```

## Example Usage
Once you make the `ldLookupx` executable, you can use it like this:

```
>>> head data.ld 
CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2
1 11008 1:11008:C:G 0.0884692 1 11012 1:11012:C:G 0.0884692 1
1 46285 1:46285:ATAT:A 0.000994036 1 66461 1:66461:T:A 0.000994036 1
1 46285 1:46285:ATAT:A 0.000994036 1 81590 rs202072409:81590:AC:A 0.000994036 1
1 49554 1:49554:A:G 0.0636183 1 76838 1:76838:T:G 0.0616302 0.889466
>>> ./ldLookupx myTable create -r 1 -K 2 -V 6 -R 8 -d " " data.ld 
>>> ./ldLookupx myTable retrieve -k 1:11008:C:G
Key: 1:11008:C:G
Values: 1:11012:C:G 
>>> ./ldLookupx myTable retrieve -k 1:46285:ATAT:A
Key: 1:46285:ATAT:A
Values: 1:66461:T:A rs202072409:81590:AC:A 
>>> ./ldLookupx myTable retrieve -k badKey
ERROR: Nonexistent key 'badKey'
>>> ./ldLookupx myTable retrieve -k 1:49554:A:G
ERROR: Nonexistent key '1:49554:A:G'
>>> ./ldLookupx badTable retrieve -k 1:11008:C:G
ERROR: Error opening file 'badTable.dht': open call failed.
```