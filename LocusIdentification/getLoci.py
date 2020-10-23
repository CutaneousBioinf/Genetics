import pandas as pd
import numpy as np
import unittest


def getLoci(threshold, path, gap):

    ## Streams data in one line at a time, but takes SIGNIFICANTLY longer to run w/ this method
    """ positions = np.array([])
        p_vals = np.array([])
        chromosomes = np.array([])
        file = open(path)
        for line in path:
            line = line.split(' ')
            if line[0] == 'chr':
                continue
            if int(line[9]) <= threshold:
                positions.append(float(line[1]))
                p_vals.append(float(line[9]))
                chromosomes.append(float(line[0])) """

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

    # stacks positions onto p values and sorts by positions
    significantMarkers = np.vstack((chromosomes, positions, p_vals))

    # finding markers with significant gaps
    mostSignificantMarkers = []
    for i in range(1, len(significantMarkers)):
        prevPosition = significantMarkers[1][i - 1]
        position = significantMarkers[1][i]
        if abs(prevPosition - position) >= gap:
            mostSignificantMarkers.append(significantMarkers[i])

    return significantMarkers, mostSignificantMarkers


def main():
    print(getLoci(0.5, "data.txt", 500000))


if __name__ == "__main__":
    main()
