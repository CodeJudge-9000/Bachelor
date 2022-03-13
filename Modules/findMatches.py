from FingerPrints import *
from atomic_annihilator import remove_atoms
from Nearest_Neighbors import *
from ase.io import read

def find_finger_print_match_simple(index, struct, removalList):
    
    system = struct.copy()
    allNNIndexes = get_removed_indeces(sys = system, theIndex = index)


    RemovalLists = make_list_of_lists(allNNIndexes)
    RemovalLists.append([])
    fingerList = []
    for rem in RemovalLists:
        copySys = struct.copy()
        atomIndex = remove_atoms(system = copySys, atomIndex = index, atomRemoveIndex = rem, relax = False, overwriteCalc = False)
        fingerList.append(SimpleFinger(system = copySys, index = atomIndex))
    ## Check whitch situation corresponds

    locVec = []

    for i in fingerList:
        if(i == removalList):
            locVec.append(1)
        else:
            locVec.append(0)
    position = locVec.index(1)
    correspondingList = RemovalLists[position]
    
    return correspondingList


def find_finger_print_match_short_common(index, struct, removalList):
    
    system = struct.copy()
    allNNIndexes = get_removed_indeces(sys = system, theIndex = index)


    RemovalLists = make_list_of_lists(allNNIndexes)
    RemovalLists.append([])
    fingerList = []
    for rem in RemovalLists:
        copySys = struct.copy()
        atomIndex = remove_atoms(system = copySys, atomIndex = index, atomRemoveIndex = rem, relax = False, overwriteCalc = False)
        fingerList.append(ShortCommonFinger(system = copySys, index = atomIndex))
    
    ## Check whitch situation corresponds

    locVec = []

    for i in fingerList:
        if(i == removalList):
            locVec.append(1)
        else:
            locVec.append(0)
    position = locVec.index(1)
    correspondingList = RemovalLists[position]
    
    return correspondingList

def search_for_case_simple(theList, sysString):
    struct = read(sysString)
    i = 0
    while True:
        try:
            firstFoundList = find_finger_print_match_simple(i, struct, theList)
            assert struct[i].symbol == 'S'
            break
        except:
            i += 1
    return i, firstFoundList

def search_for_case_short_common(theList, sysString):
    struct = read(sysString)
    i = 0
    while True:
        try:
            firstFoundList = find_finger_print_match_short_common(i, struct, theList)
            assert struct[i].symbol == 'S'
            break
        except:
            i += 1
    return i, firstFoundList