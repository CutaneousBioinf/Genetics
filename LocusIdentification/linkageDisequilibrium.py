import pandas as pd
import numpy as np


def linkageDisequilibrium(lociPath, LDPath, LDthreshold):
    locidf = pd.read_csv(lociPath)
    LDdf = pd.read_csv(LDPath)
    loci = locidf.to_numpy()
    LD = LDdf.to_numpy()

    markers = []
    # very rough draft
    for i in range(0, len(LD)):
        if (LD[i][2]) >= LDthreshold:
            markers.append(LD[i])

    return markers
