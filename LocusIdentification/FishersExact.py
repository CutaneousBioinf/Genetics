from scipy import stats
import numpy as np


def fishersExactTest(binaryFile, alternative):
    firstBEDnumones = 0
    firstBEDnumzeros = 0
    BEDnumones = 0
    BEDnumzeros = 0
    with open(binaryFile) as f:
        for line in f:
            L = line.strip().split()
            if L[2] == 0:
                firstBEDnumones = firstBEDnumones + 1
            elif L[2] == 1:
                firstBEDnumzeros = firstBEDnumzeros + 1

            for i in range(3, len(L)):
                if L[i] == 0:
                    BEDnumones = BEDnumones + 1
                elif L[i] == 1:
                    BEDnumzeros = BEDnumzeros + 1

    contingency = np.array([[firstBEDnumones, BEDnumones], [firstBEDnumzeros, BEDnumzeros]])

    p_value, oddsratio = stats.fisher_exact(contingency, alternative)

    return p_value
