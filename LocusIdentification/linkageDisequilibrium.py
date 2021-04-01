import sys
import pandas as pd
import numpy as np


def linkageDisequilibrium(lociPath, LDPath, LDthreshold, outputfile):
    locidf = pd.read_csv(lociPath)
    LDdf = pd.read_csv(LDPath)
    lociPositions = np.array(locidf["pos"].to_list())
    LD = LDdf.to_numpy()

    markers = []

    # Please present the results with the locus marker first, the LD (surrogate) marker second, and the LD value third.

    for i in range(0, len(LD)):
        if (LD[i][2]) >= LDthreshold and (LD[i][0] in lociPositions or LD[i][1] in lociPositions):
            markers.append([lociPositions[i], LD[i][0], LD[i][2]])

    markers_df = pd.DataFrame(data=markers, columns=["locusMarker", "surrogateMarker", "LDval"])
    markers_df.to_csv(outputfile) if outputfile != " " else print(markers_df)
    sys.stderr.write("Markers file created.\n")
    return markers


def main():
    print(linkageDisequilibrium("loci.txt", "ld.txt", 0.5))


if __name__ == "__main__":
    main()
