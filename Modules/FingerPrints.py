########################################################################################################################
#
# This module contain the following fingerprints (with their formats):
#
# SimpleFinger:      (Sum(Mo), Sum(S at Mo1), Sum(S at Mo2), Sum(S at Mo3))
# LessSimpleFinger:  (aoi(Mo), Mo1(S), Mo2(S), Mo3(S), aoi(S))
# CommonFinger:      (aoi(Mo), Mo1(S) and aoi(S), Mo1(S) and not aoi(S), ...)
# ShortCommonFinger: (aoi(Mo), Mo1(S) and aoi(S), Mo2(S) and aoi(S), ...)
# LongCommonFinger:  (aoi(Mo), aoi(S), Mo1(S) and aoi(S), Mo1(S) and not aoi(S), ...)
#
########################################################################################################################

# Initialization
from Nearest_Neighbors import *
from atomic_annihilator import remove_atoms

def SimpleFinger(system, index):
    """
    Fingerprint in format (Sum(Mo), Sum(S at Mo1), Sum(S at Mo2), Sum(S at Mo3), sum(S at Mo4)). If there are no Mo2, Mo3 or Mo4, then the number of neighboring S will be set to 0.
    
    Input:
     - system: Atoms ase system
     - index: which atom to calculate fingerprint for
     
    Output:
     - fingerprint: (aoi(Mo), Mo1(S), Mo2(S), Mo3(S), Mo4(S))
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
    
    # If there are less than 4 Mo atoms, add '0' to the fingerprint for each missing atom
    while(nMo < 4):
        fingerPrint.append(0)
        nMo += 1
    
    # Getting rid of degenerate situations
    first = fingerPrint[0:2]
    rest = fingerPrint[2:6]
    rest.sort(reverse=True)
    fingerPrint = first + rest
    
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
    
    # If there are less than 4 Mo atoms, add '0' to the fingerprint for each missing atom
    while(nMo < 4):
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
    
    # If there are less than 4 Mo atoms, add '0' twice to the fingerprint for each missing atom
    while(nMo < 4):
        fingerPrint.append(0)
        fingerPrint.append(0)
        nMo += 1
    
    return fingerPrint


def ShortCommonFinger(system, index):
    """
    Fingerprint in format (aoi(Mo), aoi(S), Mo1(S) and aoi(S), Mo2(S) and aoi(S), ...). If there are no Mo2 or Mo3, then the number of neighboring S will be set to 0.
    
    Input:
     - system: Atoms ase system
     - index: which atom to calculate fingerprint for
     
    Output:
     - fingerprint: (aoi(Mo), aoi(S), Mo1(S) and aoi(S), Mo2(S) and aoi(S), ...)
    """
    fingerPrint = []
    
    # Find the number of neighboring Mo atoms
    MOs = NNs_Sphere(system, index, atom_type = "Mo")
    nMo = MOs[0]
    fingerPrint.append(nMo)
    
    # Add the amount of S NNs for our aoi atom to the fingerprint
    fingerPrint.append(NNs_Sphere(system, index, atom_type = "S")[0])
    
    # Find the aoi(S) neighbor ids
    aoiNeighbors = np.array([p["index"] for p in NNs_Sphere(system, index, atom_type = "S")[2]])
    
    # Now, for each Mo atom, find the amount of S atoms that overlap with the aoi S neighbors and add to the fingerprint
    for e in MOs[2]:
        # First find all the nearest neighbor id's
        tempMoNeighbors = np.array([p["index"] for p in NNs_Sphere(system, e["index"], atom_type = "S")[2]])
        
        # Then find how many overlap
        nCommon = len(np.intersect1d(tempMoNeighbors, aoiNeighbors))
        
        # Now add this to the fingerprint
        fingerPrint.append(nCommon)
    
    # If there are less than 4 Mo atoms, add '0' to the fingerprint for each missing atom
    while(nMo < 4):
        fingerPrint.append(0)
        nMo += 1
    
    # Getting rid of degenerate situations
    first = fingerPrint[0:1]
    rest = fingerPrint[1:5]
    rest.sort(reverse=True)
    fingerPrint = first + rest
    
    return fingerPrint


def LongCommonFinger(system, index):
    """
    Fingerprint in format (aoi(Mo), aoi(S), Mo1(S) and aoi(S), Mo1(S) and not aoi(S), ...). If there are no Mo2 or Mo3, then the number of neighboring S will be set to 0.
    
    Input:
     - system: Atoms ase system
     - index: which atom to calculate fingerprint for
     
    Output:
     - fingerprint: (aoi(Mo), aoi(S), Mo1(S) and aoi(S), Mo1(S) and not aoi(S), ...)
    """
    fingerPrint = []
    
    # Find the number of neighboring Mo atoms
    MOs = NNs_Sphere(system, index, atom_type = "Mo")
    nMo = MOs[0]
    fingerPrint.append(nMo)
    
    # Add the amount of S NNs for our aoi atom to the fingerprint
    fingerPrint.append(NNs_Sphere(system, index, atom_type = "S")[0])
    
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
    
    # If there are less than 4 Mo atoms, add '0' twice to the fingerprint for each missing atom
    while(nMo < 4):
        fingerPrint.append(0)
        fingerPrint.append(0)
        nMo += 1
    
    return fingerPrint

def get_all_finger_prints_3_SCF(system, theIndex):
    """
    Function to get all unique SCF (Short common fingerprint) combinations for a three removed atom
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed

    Output:
    uniqueListFinger [list] -> The corresponding unique fingerprints encountered around the ID
    uniqueListComb [list] -> All the unique list combinations
    """

    # Collecting all fingerprints and combolist which exists
    lists = get_all_combs_3(system, theIndex)

    # Collecting unique combinations to spare out some posibilities      
    fingerList = []
    for i in lists:
        sysCopy = system.copy()
        idNew = remove_atoms(sysCopy, theIndex, i, relax = False, overwriteCalc = False)
        fingerPrint = ShortCommonFinger(sysCopy, idNew)
        fingerList.append(fingerPrint)
    return fingerList, lists

def get_all_finger_prints_2_SCF(system, theIndex):
    """
    Function to get all unique SCF (Short common fingerprint) combinations for a two removed atom
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed

    Output:
    uniqueListFinger [list] -> The corresponding unique fingerprints encountered around the ID
    uniqueListComb [list] -> All the unique list combinations
    """

    # Collecting all fingerprints and combolist which exists
    lists = get_all_combs_2(system, theIndex)

    # Collecting unique combinations to spare out some posibilities    
    fingerList = []
    for i in lists:
        sysCopy = system.copy()
        idNew = remove_atoms(sysCopy, theIndex, i, relax = False, overwriteCalc = False)
        fingerPrint = ShortCommonFinger(sysCopy, idNew)
        fingerList.append(fingerPrint)
    return fingerList, lists

def get_all_finger_prints_1_SCF(system, index):
    """
    Function to get all unique SCF (Short common fingerprint) combinations for a one removed atom
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed

    Output:
    uniqueListFinger [list] -> The corresponding unique fingerprints encountered around the ID
    uniqueListComb [list] -> All the unique list combinations
    """

    # Collecting all fingerprints and combolist which exists
    lists = get_all_combs_1(system, index)

    # Collecting unique combinations to spare out some posibilities
    fingerList = []
    for i in lists:
        sysCopy = system.copy()
        idNew = remove_atoms(sysCopy, index, i, relax = False, overwriteCalc = False)
        fingerPrint = ShortCommonFinger(sysCopy, idNew)
        fingerList.append(fingerPrint)

    return fingerList, lists


def make_unique_removals_3_SCF(system, index):
    """
    Function to get all unique SCF (Short common fingerprint) combinations for a two removed atom
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed

    Output:
    uniqueListFinger [list] -> The corresponding unique fingerprints encountered around the ID
    uniqueListComb [list] -> All the unique list combinations
    """

    # Collecting all fingerprints and combolist which exists
    fingerList, combList  = get_all_finger_prints_3_SCF(system, index)

    # Collecting unique combinations to spare out some posibilities
    uniqueListFinger = []
    uniqueListComb = []
    for i, element in enumerate(fingerList):
        if(element not in uniqueListFinger):
            uniqueListFinger.append(fingerList[i])
            uniqueListComb.append(combList[i])
    
    return uniqueListFinger, uniqueListComb

def make_unique_removals_2_SCF(system, index):
    """
    Function to get all unique SCF (Short common fingerprint) combinations for a two removed atom
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed

    Output:
    allfingers [list] -> All unique fingerprints for the unique removal lists
    allRems [list] -> All unique fingerprints for the unique removal lists
    """

    # Collecting all fingerprints and combolist which exists
    fingerList, combList = get_all_finger_prints_2_SCF(system, index)

    # Collecting unique combinations to spare out some posibilities
    uniqueListFinger = []
    uniqueListComb = []
    for i, element in enumerate(fingerList):
        if(element not in uniqueListFinger):
            uniqueListFinger.append(fingerList[i])
            uniqueListComb.append(combList[i])
            
    return uniqueListFinger, uniqueListComb

def make_unique_removals_1_SCF(system, index):
    """
    Function to get all unique SCF (Short common fingerprint) combinations for a single removed atom
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed

    Output:
    allfingers [list] -> All unique fingerprints for the unique removal lists
    allRems [list] -> All unique fingerprints for the unique removal lists
    """

    # Collecting all fingerprints and combolist which exists
    fingerList, combList = get_all_finger_prints_1_SCF(system, index)

    # Collecting unique combinations to spare out some posibilities
    uniqueListFinger = []
    uniqueListComb = []
    for i, element in enumerate(fingerList):
        if(element not in uniqueListFinger):
            uniqueListFinger.append(fingerList[i])
            uniqueListComb.append(combList[i])

    return uniqueListFinger, uniqueListComb

def get_unqiue_removal_combi(system, index):
    """Function for getting one unique combo of removal atoms in the SCF (Short common fingerprint framwork),
    which gives out one sugestions for the system and the index, and it will spit out all fingerprints it finds
    and all unique combos of either 1,2 and 3 removed atoms.
    input:
    system [ASE Object] -> the system one observes 
    index [Int] -> index of the atom which nearest neigbohrs should be removed 

    Output:
    allfingers [list] -> All unique fingerprints for the unique removal lists
    allRems [list] -> All unique fingerprints for the unique removal lists
    """
    # The making of the removals have to be in different functions because, they use different amount of for loops
    fingers1RemovedAtom, listFor1RemovedAtom = make_unique_removals_1_SCF(system, index)
    fingers2RemovedAtom, listFor2RemovedAtom = make_unique_removals_2_SCF(system, index)
    fingers3RemovedAtom, listFor3RemovedAtom = make_unique_removals_3_SCF(system, index)
    
    allFingers = fingers1RemovedAtom + fingers2RemovedAtom + fingers3RemovedAtom
    allRems = listFor1RemovedAtom + listFor2RemovedAtom + listFor3RemovedAtom

    return allFingers, allRems