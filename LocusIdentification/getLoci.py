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
@click.option('--gap', default=50000)
@click.option('--outputfile', default="loci.csv")
@click.option('--outputmarkers', default=1, help='1 to create a file of markers, and 0 to only create a '
                                                 'file for loci')
@click.option('--outputmarkersfile', default="markers.csv")
def getLoci(threshold, path, gap, chromosome, outputfile, outputmarkers, outputmarkersfile):
    # process file into csv dataframe
    df = pd.read_csv(path, delimiter='	')
    df = df[pd.to_numeric(df['p_value'], errors='coerce').notnull()]
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
    significantMarkers = significantMarkers[significantMarkers[:, 0].argsort()]
    # sort by chromosome
    significantMarkers = significantMarkers[significantMarkers[:, 1].argsort()]

    # filtering by p value
    significantMarkers = significantMarkers[np.any(significantMarkers <= threshold, axis=1), :]

    # filtering by chromosome
    if chromosome != 0:
        significantMarkers = significantMarkers[:, significantMarkers[1] != chromosome]

    sys.stderr.write("Data sorted and filtered.\n")

    significantLoci = []
    rows = np.shape(significantMarkers)[0]

    # gets loci from each chromosome, and puts them in array
    significantLocus = significantMarkers[0]
    currentChromosome = significantMarkers[0][2]
    # set initial p value to a large number
    p_val = 1
    for i in range(1, rows):
        prevPosition = significantMarkers[i - 1][0]
        position = significantMarkers[i][0]

        if significantMarkers[i][1] != currentChromosome or position - prevPosition > gap:
            currentChromosome = significantMarkers[i][1]
            p_val = 1
            significantLoci.append(significantLocus)

        if significantMarkers[i][2] <= p_val:
            significantLocus = significantMarkers[i, :]
            p_val = significantMarkers[i][2]

    sys.stderr.write("Loci identified.\n")

    if(outputmarkers):
        markers_df = pd.DataFrame(data=significantMarkers, columns=["pos", "chr", "p_value"])
        markers_df.to_csv(outputmarkersfile) if outputmarkersfile != " " else print(markers_df)
        sys.stderr.write("Markers file created.\n")

    final_df = pd.DataFrame(data=significantLoci, columns=["pos", "chr", "p_value"])
    final_df.to_csv(outputfile) if outputfile != " " else print(final_df)
    sys.stderr.write("Loci file created.\n")
    return significantLoci


def testcases():
    # one chromosome
    significantLoci = getLoci(0.5, "oneChr.txt", 50000, 6, " ")
    assert (significantLoci[0][1] == significantLoci[-1][1])

    # one of the chromosomes only contains one marker
    significantLoci = getLoci(0.5, "oneMarker.txt", 50000, 0, " ")
    assert (len(significantLoci) == 1)

    # non numeric p-val
    significantLoci = getLoci(0.5, "nonNumericP.txt", 50000, 0, " ")
    assert (type(significantLoci[0][1]) == int)


def main():
    significantLoci = getLoci()
    testcases()


if __name__ == "__main__":
    main()
