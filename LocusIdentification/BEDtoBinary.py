import numpy as np


# filenames is a list of the BED documents to compare
# markersFile is the list of markers
def convert(markersFile, filenames):
    # initial variables
    chromosomes = []
    positions = []
    binaryData = []
    headers = ["chr", "position"]

    # Read through markers file, logging the chromosomes and the position
    with open(markersFile) as f:
        for line in f:
            L = line.strip().split()
            chromosomes.append(L[0])
            positions.append(L[1])

    # Read through each remaining BED file
    for file in filenames:
        with open(file) as f:

            # Add filename to header.
            headers.append(file)
            for i in range(0, len(chromosomes)):
                found = False
                for line in f:
                    L = line.strip().split()
                    binary = []
                    # If the chromosomes being currently compared are the same, continue
                    if int(L[0]) == chromosomes[i]:
                        # If the positions of the chromosomes are overlapping, then append "1".
                        if isOverlapping(positions[i], [L[1], L[2]]):
                            binary.append(1)
                            found = True
                            break
                    else:
                        if not found:
                            binary.append(0)

            # Add binary data from this BED file to the "master" matrix.
            binaryData.append(binary)

    finalMatrix = np.vstack((chromosomes, positions, binaryData))
    np.transpose(finalMatrix)
    finalMatrix = np.vstack((headers, finalMatrix))
    return finalMatrix


def isOverlapping(position, positionIntervalArray):
    return (position > positionIntervalArray[0]) and (positionIntervalArray[1] > position)
