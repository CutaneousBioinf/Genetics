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
def getLoci(threshold, path, gap, chromosome):
    # process file into csv dataframe
    df = pd.read_csv(path, delimiter='	')
    df.sort_values(by=['pos'], inplace=True)

    # isolate p values, chr, and pos columns
    chromosomes = np.array(df["chr"].to_list())
    p_vals = np.array(df["p_value"].to_list())
    positions = np.array(df["pos"].to_list())

    # delete dataframe variable, python garbage collector deletes it from memory
    del df

    # compare to threshold
    positions = positions[p_vals <= threshold]
    chromosomes = chromosomes[p_vals <= threshold]
    p_vals = p_vals[p_vals <= threshold]

    # filter by inputted chromosome
    if chromosome != 0:
        positions = positions[chromosomes == chromosome]
        p_vals = p_vals[chromosomes == chromosome]
        chromosomes = chromosomes[chromosomes == chromosome]

    # stacks positions onto p values and sorts by positions, then by chromosomes
    significantMarkers = np.vstack((positions, chromosomes, p_vals))
    significantMarkers = significantMarkers[:, significantMarkers[0].argsort()]
    significantMarkers = significantMarkers[:, significantMarkers[1].argsort()]

    # finds the locus if a specific chromosome is defined
    significantLocus = []
    if chromosome != 0:
        currentSignificantMarker = [0, 0, 0]
        for i in range(1, len(significantMarkers)):
            prevPosition = significantMarkers[0][i - 1]
            position = significantMarkers[0][i]
            if abs(prevPosition - position) >= gap and significantMarkers[2][i] > currentSignificantMarker[2]: # p
                # value inequality
                currentSignificantMarker = significantMarkers[i]

        significantLocus.append(currentSignificantMarker)

    # gets loci from each chromosome, and puts them in array
    else:
        for j in range(1, 24):
            significantLocus = [0, 0, 0]
            for i in range(1, len(significantMarkers)):
                prevPosition = significantMarkers[0][i - 1]
                position = significantMarkers[0][i]
                if abs(prevPosition - position) >= gap and significantMarkers[2][i] > significantLocus[2] and \
                        significantMarkers[1][i] == j: # p value inequality and making sure that the marker is a
                    # chromosome that we want
                    currentSignificantMarker = significantMarkers[i]
            significantLocus.append(currentSignificantMarker)

    return significantMarkers, significantLocus


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
