import pandas as pd
import numpy as np
import click
import unittest

'''notes: This process takes a max of about 4500MB of memory while running. (for a 20mil line database)
Sorting while keeping the data in the dataframe takes longer than the current process, due to the xtra looping necessary
to remove rows with P-values under the threshold.'''


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

    # delete dataframe from memory
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

    # stacks positions onto p values and sorts by positions
    significantMarkers = np.vstack((chromosomes, positions, p_vals))

    # finding markers with significant gaps
    mostSignificantMarkers = []
    for i in range(1, len(significantMarkers)):
        prevPosition = significantMarkers[1][i - 1]
        position = significantMarkers[1][i]
        if abs(prevPosition - position) >= gap and significantMarkers[0][i - 1] == significantMarkers[0][
            i]:  # making sure 2 datapoints from different chromosomes dont get flagged
            mostSignificantMarkers.append(significantMarkers[i])

    return significantMarkers, mostSignificantMarkers


# unit tests
class Test(unittest.TestCase):

    # significant markers and most significant markers should be the same with a gap size of 0
    def testMostSig(self):
        self.assertEqual(getLoci()[0] == getLoci()[1])

    # assuming an input of 2 to the --chromosome part of CLI
    def testChromosomeIsolation(self):
        self.assertEqual(getLoci()[0][0][0] == 2)

    # come up with more tests


def main():
    print(getLoci())
    unittest.main()


if __name__ == "__main__":
    main()
