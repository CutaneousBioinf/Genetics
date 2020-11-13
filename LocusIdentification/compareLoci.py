import numpy as np
import pandas as pd

# !!! Currently still a draft
def compareLoci(path1, path2):
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2)

    loci1 = df1.to_numpy()
    loci2 = df2.to_numpy()

    similarloci = []
    for i in range(0,len(loci1)):
        for j in range(0,len(loci2)):
            if abs(loci1[i][0] - loci2[j][0]) >= 500:
                similarloci.append(loci2)

    return similarloci