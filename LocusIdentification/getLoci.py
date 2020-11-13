import pandas as pd
import numpy as np
import click

'''As of 11/6/2020: This process takes a max of about 5000MB of memory while running. (for a 20mil line database) '''


# command line interface
@click.command()
@click.option('--chromosome', default=0,
              help='chromosome to analyze, enter 0 to get all chromosome data. 0 by default')
@click.option('--path', default="data.txt", help='.txt file to analyze')
@click.option('--threshold', default=0.5, help='p_value threshold')
@click.option('--gap', default=500000)
@click.option('--totalchromosomes', default=23)
def getLoci(threshold, path, gap, chromosome, totalchromosomes):
    # process file into csv dataframe
    df = pd.read_csv(path, delimiter='	')

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

    # finds the loci if a specific chromosome is defined
    significantLoci = []
    rows = np.shape(significantMarkers)[0]
    if chromosome != 0:
        currentSignificantMarker = [0, 0, 0]
        for i in range(1, rows):
            prevPosition = significantMarkers[i - 1][0]
            position = significantMarkers[i][0]
            if abs(prevPosition - position) >= gap and significantMarkers[i][2] < currentSignificantMarker[2]:  # p
                # value inequality
                currentSignificantMarker = significantMarkers[i, :]

        significantLoci.append(currentSignificantMarker)

    # gets loci from each chromosome, and puts them in array
    # TODO: Figure out more efficient algorithm that only uses one loop
    else:
        for j in range(1, totalchromosomes + 1):
            significantLocus = [0, 0, 0]
            for i in range(1, rows):
                prevPosition = significantMarkers[i - 1][0]
                position = significantMarkers[i][0]
                if abs(prevPosition - position) >= gap and significantMarkers[i][2] < significantLocus[2] and \
                        significantMarkers[i][1] == j:  # p value inequality and making sure that the marker is a
                    # chromosome that we want
                    significantLocus = significantMarkers[i, :]
            significantLoci.append(significantLocus)

    final_df = pd.DataFrame(data=significantLoci, columns=["pos", "chr", "p_value"])
    final_df.to_csv("loci.csv")
    return significantMarkers, significantLoci


def main():
    significantMarkers, significantLocus = getLoci()

    # makes sure that there are at most 23 loci, one for each chromosome
    assert len(significantLocus) <= 23

    # There can never be more loci than markers
    assert len(significantMarkers) >= len(significantLocus)

    # true if positions are sorted
    assert significantMarkers[0][0] < significantMarkers[0][1]

    print(significantMarkers, significantLocus)


if __name__ == "__main__":
    main()
