from FingerPrints import *
from atomic_annihilator import remove_atoms
from Nearest_Neighbors import *

def find_finger_print_match_simple(index, struct, removalList):
    
    system = struct.copy()
    allNNIndexes = get_removed_indeces(sys = system, theIndex = index)


    RemovalLists = make_list_of_lists(allNNIndexes)
    RemovalLists.append([])
    fingerList = []
    for rem in RemovalLists:
        copySys = struct.copy()
        atomIndex = remove_atoms(system = copySys, atomIndex = index, atomRemoveIndex = rem, relax = False)
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
        atomIndex = remove_atoms(system = copySys, atomIndex = index, atomRemoveIndex = rem, relax = False)
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
