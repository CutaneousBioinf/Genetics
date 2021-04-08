import numpy as np


# filenames is a list of the BED documents to compare
# lociFile is the list of loci
def convert(markersFile, BEDfilenames):
    # initial variables
    chromosomes = []
    lociPositions = []
    positions = []
    binaryData = []
    headers = ["chr", "pos"]

    with open(markersFile) as f:
        for line in f:
            L = line.strip().split()
            chromosomes.append(L[0])
            lociPositions.append(L[3])
            positions.append(L[1])  # (chr, pos, lociPos)

    loci = np.unique(lociPositions).sort()

    # Read through each remaining BED file
    for file in BEDfilenames:
        with open(file) as f:
            # Add filename to header.
            headers.append(file)
            for j in range(loci):
                for i in range(0, len(positions)):
                    found = False
                    for line in f:
                        L = line.strip().split()
                        binary = []
                        # If the chromosomes being currently compared are the same, continue
                        if int(L[0]) == chromosomes[i] and loci[j] == lociPositions[i]:
                            if isOverlapping(lociPositions[i], [L[1], L[2]]) or isOverlapping(positions[i], [L[1], L[2]]):
                                binary.append(1)
                                found = True
                                break;
                if not found:
                    binary.append(0)

            # Add binary data from this BED file to the "master" matrix.
            binaryData.append(binary)

    finalMatrix = np.vstack((chromosomes, loci, binaryData))
    np.transpose(finalMatrix)
    finalMatrix = np.vstack((headers, finalMatrix))
    return finalMatrix


def isOverlapping(position, positionIntervalArray):
    return (position > positionIntervalArray[0]) and (positionIntervalArray[1] > position)
