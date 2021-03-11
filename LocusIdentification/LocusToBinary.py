import numpy as np


# Things to ask:
# Is the data going in as files or as lists/dataframes?
# How do you know a marker is in a locus? Is it based on position or chr?
# Should the significant markers be saved in getLoci.py?
def convert(markersFile, lociList):
    chromosomes = []
    positions = []
    binaryData = []
    headers = ["chr", "position"]

    # Go through each loci, logging the chromosomes and the position
    with open(lociList) as f:
        for line in f:
            L = line.strip().split()
            chromosomes.append(L[0])
            positions.append(L[1])

    # Read through each remaining BED file
    with open(markersFile) as f:
        for i in range(0, len(positions)):
            # From here, mostly the same as BED to binary
            found = False
            for line in f:
                L = line.strip().split()
                binary = []
                # add the position of the marker to the headers
                headers.append(L[1])
                # If the chromosomes being currently compared are the same, continue
                if int(L[0]) == chromosomes[i]:
                    # If the positions of the markers and the loci are overlapping, then append "1".
                    if isOverlapping(positions[i], L[1]):
                        binary.append(1)
                        found = True
                        break
            if not found:
                binary.append(0)

            # Add binary data from this marker to the "master" matrix.
            binaryData.append(binary)

    finalMatrix = np.vstack((chromosomes, positions, binaryData))
    np.transpose(finalMatrix)
    finalMatrix = np.vstack((headers, finalMatrix))
    return finalMatrix


# not sure about this
def isOverlapping(position, markerPosition):
    return position is markerPosition
