import sys
import pandas as pd
import numpy as np


def linkageDisequilibrium(lociPath, LDPath, LDthreshold, outputfile):
    locidf = pd.read_csv(lociPath)
    LDdf = pd.read_csv(LDPath)
    markerPositions = np.array(locidf["pos"].to_list())
    markerChr = np.array(locidf["chr"].to_list())
    markerPVal = np.array(locidf["p_value"].to_list())
    lociPos = np.array(locidf["correspondingLociPos"].to_list())
    LD = LDdf.to_numpy()

    markers = []

    for i in range(0, len(LD)):
        if (LD[i][2]) >= LDthreshold:
            for j in range(0, len(markerPositions)):
                if LD[i][0] == j or LD[i][1] == j:
                    markers.append([markerPositions[j], markerChr[j], markerPVal[j], lociPos[j], LD[i][2]])

    markers_df = pd.DataFrame(data=markers, columns=["pos", "chr", "p_value", "correspondingLociPos", "LDval"])
    markers_df.to_csv(outputfile) if outputfile != " " else print(markers_df)
    sys.stderr.write("Markers file created.\n")
    return markers


def main():
    print(linkageDisequilibrium("loci.txt", "ld.txt", 0.5))


if __name__ == "__main__":
    main()
