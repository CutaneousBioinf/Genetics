# rsLookup

The rsLookup library allows for the searching of the snp150Common dataset for more efficient manipulation of genomic markers.

## Installation

Simply clone the repository and run the make command, which will generate a self-contained rsLookup binary. The diskhash library is included and will be statically linked with the rsLookup binary.

## Commands to prepare and reference hash tables

To create the rsID -> chromosome/position/allele map:
`rsLookup -C -k rsid -s [path_to_snp150Common.txt] -t [rsid_table_name.dht]`
To create the chromosome/position/allele -> rsID map:
`rsLookup -C -k cpa -s [path_to_snp150Common.txt] -t [cpa_table_name.dht]`
To return the rsID of a marker given the chromosome/position/allele:
`rsLookup -R -k rsid -t [rsid_table_name.dht] -r rsNNNNNNNNN`
To return the rsID of a marker given the chromosome/position/allele sequence:
`rsLookup -R -k cpa -t [rsid_table_name.dht] -c chrXX -p NNNNN -a [allele_sequence]`
Note: Positions are indexed from 0 when returned but indexed from 1 when supplied by the user.

## Other options
- The `-e` option outputs an error log to standard error, which will include information about markers that do not match user input during lookups and also markers that are not included or recognized in the hash table creation process. Alternatively, the user may supply `-l` followed by a filename to send the same error log to a file.
- The `-f` option allows for input from a file. For rsIDs, they should be one each on their own line. Each line of a reverse lookup should have a tab-separated string of the chromosome, position, and allele on each line.

