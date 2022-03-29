from FingerPrints import *
from atomic_annihilator import remove_atoms
from Nearest_Neighbors import *
from ase.io import read
"""
A package for finding matches to missing fingerprints encountered when doing the KMC - Simulations.
"""

def find_finger_print_match_simple(index, struct, fingerPrint):
    """
    A function that tries to find a combination of indexes to remove to gain the simple fingerprint one searches for
    
    Input:
    ------
    - `index` [Int] -> the atom one wants to move
    - `struct` [ASE] -> The system we consider
    - `fingerPrint` [list] -> A fingerprint list of the type (Simple Finger)

    Output:
    -------
     - `correspondingList` [list] -> The list with the indexes we should remove
    The removalist, that gives the fingerprint
    """
    
    # Getting the system, and all possible indexes to remove
    system = struct.copy()
    allNNIndexes = get_removed_indeces(system, index)

    # Using the make_list_of_lists function to go from index list to a list of listst to remove and appending with the case where nothing got removed.
    RemovalLists = make_list_of_lists(allNNIndexes)
    RemovalLists.append([])

    # Removing all the atoms in the removal list, to get their fingerprint and add them to a list
    fingerList = []
    for rem in RemovalLists:
        copySys = struct.copy()
        atomIndex = remove_atoms(system = copySys, atomIndex = index, atomRemoveIndex = rem, relax = False, overwriteCalc = False)
        fingerList.append(SimpleFinger(system = copySys, index = atomIndex))

    # Make a logical list of type [0,0,0,0,0,0,1....] To then find the first occuring 1 which gives the first removalist that gives the fingerprint one looks for 
    locVec = []
    for i in fingerList:
        if(i == fingerPrint):
            locVec.append(1)
        else:
            locVec.append(0)
    # This line is where the function might give and error which is fine, sometimes we are not guaranteed that our fingerprint can exist for the index
    position = locVec.index(1)
    correspondingList = RemovalLists[position]
    
    return correspondingList


def find_finger_print_match_short_common(index, struct, fingerPrint):
    """
    A function that tries to find a combination of indexes to remove to gain the short common fingerprint one searches for
   
    Input:
    ------
     - `index` [Int] -> the atom one wants to move
     - `struct` [ASE] -> The system we consider
     - `fingerPrint` [list] -> A fingerprint list of the type (short common)

    Output:
    -------
     - `correspondingList` [list] -> The list with the indexes we should remove
    The removalist, that gives the fingerprint
    """

    # Getting the system, and all possible indexes to remove       
    system = struct.copy()

    
    # Using the make_list_of_lists function to go from index list to a list of listst to remove and appending with the case where nothing got removed.
    #### Old functionality
    ####allNNIndexes = get_removed_indeces(sys = system, theIndex = index)
    ####RemovalLists = make_list_of_lists(allNNIndexes)
    ####RemovalLists.append([])
    
    #### New functionality
    RemovalLists = get_all_combs_tot(system, index, True) # True for degenerate is also searched within
    RemovalLists.append([])

    # Removing all the atoms in the removal list, to get their fingerprint and add them to a list
    fingerList = []
    for rem in RemovalLists:
        copySys = struct.copy()
        atomIndex = remove_atoms(system = copySys, atomIndex = index, atomRemoveIndex = rem, relax = False, overwriteCalc = False)
        fingerList.append(ShortCommonFinger(system = copySys, index = atomIndex))

    # Make a logical list of type [0,0,0,0,0,0,1....] To then find the first occuring 1 which gives the first removalist that gives the fingerprint one looks for 
    locVec = []
    for i in fingerList:
        if(i == fingerPrint):
            locVec.append(1)
        else:
            locVec.append(0)
    # This line is where the function might give and error which is fine, sometimes we are not guaranteed that our fingerprint can exist for the index
    position = locVec.index(1)
    correspondingList = RemovalLists[position]
    
    return correspondingList



def search_for_case_simple(fingerPrintToFind, sysString):
    """
    A Function to search the ASE object for a simple finger fingerprint, to find the first occurence
    
    Input:
    ------
     - `fingerPrintToFind` [list] -> The simple fingerprint we want find
     - `sysString` [Str] -> String of the trajectory file we want to find the fingerprint in

    Output:
    -------
     - `atomIndex` [Int] -> Index of the atom we found with the fingerPrint
     - `firstFoundRemovalList` [list] -> The list of the atom indexes we want to remove to gain the fingerprint we sought

    """
    system = read(sysString)

    # Making a loop that inserts the list and tries to use the function in all indexes, on the system if we get an error because we did not find it, we will \
    # instead try on the next index in the structure until we hit something is in aggreement with our list or we have run out of indeces
    atomIndex = 0
    while True:
        if(atomIndex > len(system) + 1 ):
            raise UserWarning("You have done something wrong!")
            break
        try:
            firstFoundRemovalList = find_finger_print_match_simple(atomIndex, system, fingerPrintToFind)
            assert system[atomIndex].symbol == 'S'
            break
        except:
            atomIndex += 1

    return atomIndex, firstFoundRemovalList


def search_for_case_short_common(fingerPrintToFind, sysString, startFrom = 0):
    """
    A Function to search the ASE object for a short common finger fingerprint, to find the first occurence
    
    Input:
    ------
     - `fingerPrintToFind` [list] -> The short common fingerprint we want find
     - `sysString` [Str] -> String of the trajectory file we want to find the fingerprint in
     - `startFromt` [int] -> the atom we want to start from

    Output:
    -------
     - `atomIndex` [Int] -> Index of the atom we found with the fingerPrint
     - `firstFoundRemovalList` [list] -> The list of the atom indexes we want to remove to gain the fingerprint we sought

    """
    system = read(sysString)

    # Making a loop that inserts the list and tries to use the function in all indexes, on the system if we get an error because we did not find it, we will \
    # instead try on the next index in the structure until we hit something is in aggreement with our list or we have run out of indeces
    atomIndex = startFrom
    while True:
        if(atomIndex > len(system) + 1 ):
            raise UserWarning("You have done something wrong!")
            break
        try:
            firstFoundList = find_finger_print_match_short_common(atomIndex, system, fingerPrintToFind)
            assert system[atomIndex].symbol == 'S'
            break
        except:
            atomIndex += 1
    return atomIndex, firstFoundList

def extensive_search_short_common(fingerPrintToFind, sysString, numberOfCases):
    """ extensive_search_short_common - Is used to find the amount of cases we want in our fingerprint
     
     Inputs:
     -------
      - `fingerPrintToFind` [list] -> the list of the short common finger we want to search for
      - `sysString` [str] -> The string of the system we investigate
      - `numberOfCases` [int] -> The number of different rem lists we want to achieve a fingerprint
     
     Output:
     -------
      - `casesFound` [list] -> A list with the tuples that the search_for_case_short_common returns.
          in format of [ (int, remlist), (int, remlist)......] and so on
    
    """
    casesFound = []
    index = 0 
    
    # We look for a case, and when we find it, we start from the next atom in the sequence and find it
    for case in range(numberOfCases):
        index, remList = search_for_case_short_common(fingerPrintToFind, sysString, index)
        casesFound.append((index, remList))
        index = index + 1
    return casesFound


def get_all_requests_simple(searchList, sysString):
    """
    A function to get all the requested fingerprint in the (simple finger) format, and give out a list of indexes with removalists, that contains that.
    
    Input:
    ------
     - `searchList` [list] -> A list of the fingerprints that we want to find
     - `sysString` [Str] -> String of the trajectory file we want to find the fingerprint in

    Output:
    -------
     - `elementList` [list] -> A list with the elements of type [index, removList], which can be unpacked to so one can quickly find the fingerprints.
    """

    elementList = []
    for i in range(0,len(searchList)):
        theElement = searchList[i]
        index, RemovList = search_for_case_simple(theElement, sysString)
        element = [index, RemovList]
        elementList.append(element)
        
    return elementList