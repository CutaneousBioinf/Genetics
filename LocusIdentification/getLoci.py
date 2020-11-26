import pandas as pd
import numpy as np
import click
import sys

'''As of 11/6/2020: This process takes a max of about 5000MB of memory while running. (for a 20mil line database) '''


# command line interface
@click.command()
@click.option('--chromosome', default=0,
              help='chromosome to analyze, enter 0 to get all chromosome data. 0 by default')
@click.option('--path', default="data.txt", help='.txt file to analyze')
@click.option('--threshold', default=0.5, help='p_value threshold')
@click.option('--gap', default=500000)
@click.option('--outputfile', default=" ")
def getLoci(threshold, path, gap, chromosome, outputfile):
    # process file into csv dataframe
    df = pd.read_csv(path, delimiter='	')
    sys.stderr.write("File parsed.\n")

    # isolate p values, chr, and pos columns
    chromosomes = np.array(df["chr"].to_list())
    p_vals = np.array(df["p_value"].to_list())
    positions = np.array(df["pos"].to_list())

    # delete dataframe variable, python garbage collector deletes it from memory
    del df

    # stacks positions onto p values and sorts by positions, then by chromosomes
    significantMarkers = np.vstack((positions, chromosomes, p_vals))
    # transposes data, so that it is easier to sort, filter, and put into csv.
    significantMarkers = np.transpose(significantMarkers)
    # positions stay sorted in place
    significantMarkers = significantMarkers[significantMarkers[0].argsort()]
    significantMarkers = significantMarkers[significantMarkers[1].argsort()]

    # filtering by p value
    significantMarkers = significantMarkers[:, significantMarkers[2] <= threshold]

    # filtering by chromosome
    if chromosome != 0:
        significantMarkers = significantMarkers[:, significantMarkers[1] != chromosome]

    sys.stderr.write("Data sorted and filtered.\n")

    significantLoci = []
    rows = np.shape(significantMarkers)[0]
    # gets loci from each chromosome, and puts them in array
    significantLocus = [0, 0, 0]
    currentChromosome = 1
    # set initial p value to a large number
    p_val = 1
    for i in range(1, rows):
        prevPosition = significantMarkers[i - 1][0]
        position = significantMarkers[i][0]

        if (abs(prevPosition - position) >= gap and significantMarkers[i][1] == currentChromosome and
            significantMarkers[i][2] < p_val) or significantMarkers[i][1] != currentChromosome:
            significantLocus = significantMarkers[i, :]

        if significantMarkers[i][1] != currentChromosome:
            currentChromosome += significantMarkers[i][1]
            p_val = 1

        significantLoci.append(significantLocus)

    sys.stderr.write("Loci identified.\n")

    final_df = pd.DataFrame(data=significantLoci, columns=["pos", "chr", "p_value"])
    final_df.to_csv(outputfile) if outputfile != " " else print(final_df)
    sys.stderr.write("Output file created.\n")
    return significantMarkers, significantLoci


def testcases(significantMarkers, significantLocus):
    # one chromosome
    assert(significantMarkers[0][1] == significantMarkers[-1][1])
    # one of the chromosomes only contains one marker
    assert (len(significantMarkers) == 1)
    # non numeric chromosomes
    # p value is N/A


def main():
    significantMarkers, significantLocus = getLoci()
    testcases(significantMarkers, significantLocus)


if __name__ == "__main__":
    main()
