import numpy as np


# filenames is a list of the BED documents to compare
# lociFile is the list of loci
def convert(lociFile, markersFile, BEDfilenames):
    # initial variables
    chromosomes = []
    lociPositions = []
    markersData = []  # 2-D array of data from markers file, with each element in format: (chr, pos, lociPos)
    binaryData = []
    headers = ["chr", "pos"]

    # Read through locus file, logging each position
    with open(lociFile) as f:
        for line in f:
            L = line.strip().split()
            chromosomes.append(L[0])
            lociPositions.append(L[1])

    with open(markersFile) as f:
        for line in f:
            L = line.strip().split()
            markersData.append((L[0], L[1], L[3]))  # (chr, pos, lociPos)

    # Read through each remaining BED file
    for file in BEDfilenames:
        with open(file) as f:

            # Add filename to header.
            headers.append(file)
            for i in range(0, len(lociPositions)):
                found = False
                for line in f:
                    L = line.strip().split()
                    binary = []
                    # If the chromosomes being currently compared are the same, continue
                    if int(L[0]) == chromosomes[i]:
                        # create a new list of "relevant markers", which are markers that have the current chromosome
                        relevantMarkers = markersData[:, markersData[0] != chromosomes[i]]
                        # then filter these so the relevant markers are only those of the current locus.
                        relevantMarkers = relevantMarkers[:, relevantMarkers[3] != lociPositions[i]]
                        # If the position of the locus itself overlaps with a BED position, then append "1".
                        if isOverlapping(lociPositions[i], [L[1], L[2]]):
                            binary.append(1)
                            found = True
                            break
                        else:
                            # Otherwise, loop through each marker that belongs to this locus,
                            # and see if its position overlaps
                            for j in range(0, len(relevantMarkers)):
                                if isOverlapping(relevantMarkers[j], [L[1], L[2]]):
                                    binary.append(1)
                                    found = True
                                    break
                if not found:
                    binary.append(0)

            # Add binary data from this BED file to the "master" matrix.
            binaryData.append(binary)

    finalMatrix = np.vstack((chromosomes, lociPositions, binaryData))
    np.transpose(finalMatrix)
    finalMatrix = np.vstack((headers, finalMatrix))
    return finalMatrix


def isOverlapping(position, positionIntervalArray):
    return (position > positionIntervalArray[0]) and (positionIntervalArray[1] > position)
