# rsLookup

The rsLookup library allows for the searching of the snp150Common dataset for more efficient manipulation of genomic markers.

## Installation

Simply clone the repository and run the make command, which will generate a self-contained rsLookup binary. The diskhash library is included and will be statically linked with the rsLookup binary. rsLookup depends on the boost-iostreams and boost-filesystem libraries.

## Creating hash tables

To create the rsID -> chromosome/position/allele map:
`rsLookup -C -k rsid -s [path_to_snp150Common.txt[.gz]] -t [rsid_table_name.dht]`
To create the chromosome/position/allele -> rsID map:
`rsLookup -C -k cpa -s [path_to_snp150Common.txt[.gz]] -t [cpa_table_name.dht]`
Alternatively, one may create both hash tables by specifying the -d option as a directory rather than the -t option.
Shortform version:
`rsLookup -C [path_to_snp150Common.txt[.gz]] [target_directory_path]`

## Retrieving from hash tables

To return the rsID of a marker given the chromosome/position/allele:
`rsLookup -R -k rsid -t [rsid_table_name.dht] -r rsNNNNNNNNN`
To return the rsID of a marker given the chromosome/position/allele sequence:
`rsLookup -R -k cpa -t [cpa_table_name.dht] -c chrXX -p NNNNN -a [allele_sequence]`
Alternatively, one may use the -d option to specify the data directory rather than the -t option.
Shortform version:
`rsLookup [data_directory_path_or_appropriate_table] rsNNNNNNNNN`
`rsLookup [data_directory_path_or_appropriate_table] chrXX [position] [allele_sequence]`
Note: One may omit the chr at the beginning of the chromosome string and expect the same behavior.

## Other options
- The `-v` option outputs an error log to standard error, which will include information about markers that do not match user input during lookups and also markers that are not included or recognized in the hash table creation process. Alternatively, the user may supply `-l` followed by a filename to send the same error log to a file.
- The `-q` option suppresses default error output for the retrieval.
- The `-f` option allows for input from a file. For rsIDs, they should be one each on their own line. Each line of a reverse lookup should have a tab-separated string of the chromosome, position, and allele on each line.

