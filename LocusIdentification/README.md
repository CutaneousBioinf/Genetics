# Locus Identification tool

This tool will identify genome-wide and suggestively significant loci from GWAS summary statistics, based on their p-value and distance from each other.

For an illustration of our goal, see this Manhattan plot: https://en.wikipedia.org/wiki/Manhattan_plot#/media/File:Manhattan_Plot.png

## Code Workflow Diagram
![Flowchart of our Code System](https://user-images.githubusercontent.com/20675645/115828739-1e373800-a3dc-11eb-9d4a-a7cd90c4234e.png "Optional title")

## getLoci.py
 This script utilizes a user defined p value and position threshold to sort through a list of summary statistics and identify significant loci and markers. It will return this list of loci, or a list of markers derived from these statistics. 
### Inputs
- Path: A text (.txt) file containing GWAS summary statistics, specifically with columns "chr", "p_value", and "pos".
- Threshold: The maximum p value that a marker must have to be considered significant.
- Gap: The minimum distance two markers must have to be considered a loci.
- Chromosome: 0 if you want to analyze all chromosomes. If its another number, then it will only analyze that specific chromosome. 
- Output File: Filename to output to. Will print to the console if empty ("").
- Output Markers: 1 to create a file of markers, and 0 to only create a file of loci. 
### Outputs
If the "Output Markers" parameter is 0: A file containing a list of all markers, along with the loci they belong to. Will contain columns "pos", "chr", "p_value", and "correspondingLociPos", specifically in that order.

If the "Output Markers" parameter is 1: A file containing a list of all loci. Will contain columns "pos", "chr",and "pos", specifically in that order. 
### Code flows into:

- linkageDisequilibrium.py
- compareLoci.py
- LocusToBinary.py
- BEDToBinary.py

## linkageDisequilibrium.py
 An alternative approach to getLoci.py. Uses a user defined minimum linkage disequilibrium value to identify and return significant markers. 
### Inputs
 - Loci Path: Path to a text (.txt) file containing a list of loci, with the first three columns containing position, chromosome,and p-value in that order. 
 - LD Path: Path to a text (.txt) file containing a list of LD values, with the first three columns containing the position of marker one, position of marker 2, and the LD value, in that order. 
 - LD Threshold: The minimum LD value for a loci to be considered significant. 
 - Output File: The file to output to. If empty (""), it will print the results to the console. 
### Outputs
 A file containing a list of significant markers, with columns "pos", "chr", "p_value", "correspondingLociPos", and "LDval". 
### Code flows into:

- compareLoci.py
- LocusToBinary.py
- BEDToBinary.py

## compareLoci.py
 Compares two sets of loci, returning a list of loci that are present in both sets, have the same chromosome, and have approximately the same position. 
### Inputs
- Path 1: Path to first text (.txt) file containing a list of loci, containing columns for position and chromosome. 
- Path 2: Path to second text (.txt) file containing a list of loci, containing columns for position and chromosome.
- Gap: The position of two loci must be equal to or less than this position to be considered equal.
### Outputs
List of similar loci, with columns "pos", "chr", and "p-value".

## BEDToBinary.py
 functionality
### Inputs
### Outputs
### Code flows into:

- FishersExact.py

## LocusToBinary.py
 functionality
### Inputs
### Outputs
### Code flows into:

- FishersExact.py

## FishersExact.py
 functionality
### Inputs
### Outputs

