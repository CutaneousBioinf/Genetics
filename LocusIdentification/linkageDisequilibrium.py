import pandas as pd
import numpy as np


def linkageDisequilibrium(lociPath, LDPath, LDthreshold):
    locidf = pd.read_csv(lociPath)
    LDdf = pd.read_csv(LDPath)
    lociPositions = np.array(locidf["pos"].to_list())
    LD = LDdf.to_numpy()

    markers = []

    for i in range(0, len(LD)):
        if (LD[i][2]) >= LDthreshold and (LD[i][0] in lociPositions or LD[i][1] in lociPositions):
            markers.append(lociPositions, LD[i][0], LD[i][2])

    return markers


def main():
    print(linkageDisequilibrium("loci.txt", "ld.txt", 0.5))


if __name__ == "__main__":
    main()
