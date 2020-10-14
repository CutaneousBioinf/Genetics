import pandas as pd
import numpy as np
import math


def getLoci(threshold, path, gap):
    # process file into csv dataframe
    df = pd.read_csv(path, delimiter='	')

    # isolate p values and pos columns
    p_vals = np.array(df["p_value"].to_list())
    positions = np.array(df["pos"].to_list())

    # compare to threshold
    positions = positions[p_vals <= threshold]
    p_vals = p_vals[p_vals <= threshold]

    # stacks positions onto p values
    significantMarkers = np.vstack((positions, p_vals))

    # finding markers with significant gaps
    mostSignificantMarkers = []
    for i in range(1, len(significantMarkers)):
        prevPosition = significantMarkers[0][i-1]
        position = significantMarkers[0][i]
        if abs(prevPosition - position) >= gap:
            mostSignificantMarkers.append(significantMarkers[i])

    return significantMarkers, mostSignificantMarkers
