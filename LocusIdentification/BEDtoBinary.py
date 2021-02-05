import numpy as np


# filenames is a list of the BED documents to compare
# regulatory is the name of the BED file that contains the regulatory regions.
def convert(regulatory, filenames):
    # initial variables
    chromosomes = []
    positions = []
    binaryData = []
    headers = ["chr", "position"]

    # Read through regulatory BED file, logging the chromosome column and the position interval
    with open(regulatory) as f:
        for line in f:
            L = line.strip().split()
            chromosomes.append(L[0])
            positions.append([L[1], L[2]])

    # Read through each remaining BED file
    for file in filenames:
        with open(file) as f:
            i = 0

            # Add filename to header.
            headers.append(file)
            binary = []
            for line in f:
                L = line.strip().split()

                # If the positions of the chromosomes are overlapping, then append "1".
                # This method assumes that all BED files have the same "chr" column.
                if isOverlapping(positions[i], [L[1], L[2]]):
                    binary.append(1)
                else:
                    binary.append(0)
                i = i + 1
            # Add binary data from this BED file to the "master" matrix.
            binaryData.append(binary)

    finalMatrix = np.vstack((chromosomes, positions, binaryData))
    np.transpose(finalMatrix)
    finalMatrix = np.vstack((headers, finalMatrix))
    return finalMatrix


def isOverlapping(array1, array2):
    return (array1[1] > array2[0]) and (array2[1] > array1[0])
