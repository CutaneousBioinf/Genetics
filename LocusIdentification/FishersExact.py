from scipy import stats
import numpy as np
import pandas as pd


def fishersExactTest(binaryFile, alternative):
    numFirstFileOnes = 0
    numFirstFileZeros = 0
    numBothOnes = 0
    numBothZeros = 0
    p_values = []
    oddsratios = []
    filenames = []
    with open(binaryFile) as f:
        for line in f:
            L = line.strip().split()
            for i in range(3, len(L)):
                if L[i] == 0 and L[2] == 0:
                    numBothZeros = numBothZeros + 1
                elif L[i] == 1 and L[2] == 1:
                    numBothOnes = numBothOnes + 1
                elif L[i] == 0 and L[2] == 1:
                    numFirstFileOnes = numFirstFileOnes + 1
                elif L[i] == 1 and L[2] == 0:
                    numFirstFileZeros = numFirstFileZeros + 1
                elif L[0] == "chr":
                    filenames.append(L[i])

            contingency = np.array([[numBothOnes, numFirstFileOnes], [numFirstFileZeros, numBothZeros]])
            p_value, oddsratio = stats.fisher_exact(contingency, alternative)

            p_values.append(p_value)
            oddsratios.append(oddsratio)

    # return a dataframe with filename and p value
    finalDataset = np.vstack((p_values, oddsratios))
    final_df = pd.DataFrame(data=finalDataset, columns=filenames)
    return final_df
