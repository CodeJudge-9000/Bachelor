from Nearest_Neighbors import *


def SimpleFinger(system, index):
    """
    Fingerprint in format (Sum(Mo), Sum(S at Mo1), Sum(S at Mo2), Sum(S at Mo3)). If there are no Mo2 or Mo3, then the number of neighboring S will be set to 0.
    
    Input:
     - system: Atoms ase system
     - index: which atom to calculate fingerprint for
     
    Output:
     - fingerprint: (aoi(Mo), Mo1(S), Mo2(S), Mo3(S))
    """
    fingerPrint = []
    
    # Find the number of neighboring Mo atoms
    MOs = NNs_Sphere(system, index, atom_type = "Mo")
    nMo = MOs[0]
    fingerPrint.append(nMo)
    
    # For each neighboring Mo atom, check how many neighoring S atoms there are (and subtract one, as to not count the atom we view)
    for e in MOs[2]:
        toAppend = NNs_Sphere(system, e["index"], atom_type = "S")[0] - 1
        fingerPrint.append(toAppend)
    
    # If there are less than 3 Mo atoms, add '0' to the fingerprint for each missing atom
    while(nMo < 3):
        fingerPrint.append(0)
        nMo += 1
    
    return fingerPrint



def LessSimpleFinger(system, index):
    """
    Fingerprint in format (aoi(Mo), Mo1(S), Mo2(S), Mo3(S), aoi(S)). If there are no Mo2 or Mo3, then the number of neighboring S will be set to 0.
    
    Input:
     - system: Atoms ase system
     - index: which atom to calculate fingerprint for
     
    Output:
     - fingerprint: (aoi(Mo), Mo1(S), Mo2(S), Mo3(S), aoi(S))
    """
    fingerPrint = []
    
    # Find the number of neighboring Mo atoms
    MOs = NNs_Sphere(system, index, atom_type = "Mo")
    nMo = MOs[0]
    fingerPrint.append(nMo)
    
    # For each neighboring Mo atom, check how many neighoring S atoms there are (and subtract one, as to not count the atom we view)
    for e in MOs[2]:
        toAppend = NNs_Sphere(system, e["index"], atom_type = "S")[0] - 1
        fingerPrint.append(toAppend)
    
    # If there are less than 3 Mo atoms, add '0' to the fingerprint for each missing atom
    while(nMo < 3):
        fingerPrint.append(0)
        nMo += 1
    
    # Add the number of aoi S neighbors
    fingerPrint.append(NNs_Sphere(system, index, atom_type = "S")[0])
    
    return fingerPrint


def CommonFinger(system, index):
    """
    Fingerprint in format (aoi(Mo), Mo1(S) and aoi(S), Mo1(S) and not aoi(S), ...). If there are no Mo2 or Mo3, then the number of neighboring S will be set to 0.
    
    Input:
     - system: Atoms ase system
     - index: which atom to calculate fingerprint for
     
    Output:
     - fingerprint: (aoi(Mo), Mo1(S) and aoi(S), Mo1(S) and not aoi(S), ...)
    """
    fingerPrint = []
    
    # Find the number of neighboring Mo atoms
    MOs = NNs_Sphere(system, index, atom_type = "Mo")
    nMo = MOs[0]
    fingerPrint.append(nMo)
    
    # Find the aoi(S) neighbor ids
    aoiNeighbors = np.array([p["index"] for p in NNs_Sphere(system, index, atom_type = "S")[2]])
    
    # Now, for each Mo atom, find the amount of S atoms that overlap with the aoi S neighbors and add to the fingerprint
    for e in MOs[2]:
        # First find all the nearest neighbor id's
        tempMoNeighbors = np.array([p["index"] for p in NNs_Sphere(system, e["index"], atom_type = "S")[2]])
        
        # Then find how many overlap, and not
        nCommon = len(np.intersect1d(tempMoNeighbors, aoiNeighbors))
        nNotCommon = len(tempMoNeighbors) - 1 - nCommon # Remember to subtract the atom itself
        
        # Now add this to the fingerprint
        fingerPrint.append(nCommon)
        fingerPrint.append(nNotCommon)
    
    # If there are less than 3 Mo atoms, add '0' twice to the fingerprint for each missing atom
    while(nMo < 3):
        fingerPrint.append(0)
        fingerPrint.append(0)
        nMo += 1
    
    return fingerPrint

