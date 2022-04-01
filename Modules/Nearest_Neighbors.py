from array import array
from multiprocessing.dummy import Array
import numpy as np
from math import sqrt
from ase import Atoms

### Function to get distance between atoms ###
def atomDist(atom1, atom2):
    d1 = atom1.position[0] - atom2.position[0]
    d2 = atom1.position[1] - atom2.position[1]
    d3 = atom1.position[2] - atom2.position[2]
    dist = sqrt(d1**2 + d2**2 + d3**2)
    return dist

#
# Note for these functions: For whole system: This-> O(2*n), If specialized -> O(n)
#
# Using these functions for a system means calculating all distances twice
#

# First function, sphere
def NNs_Sphere(system, ind, atom_type = "*Insert Symbol Here*"):
    """
    This function calculates the distances from the given atom, and all other atoms in the system. It then takes the shortest distance and uses it to draw a sphere 115% of this distance around the given atom. Everything within this sphere is counted as a neighboring atom.
    
    
    Input:
     - ase atoms system [Atoms]
     - index [int]
     - atomic symbol/type [string]
    
    Output:
     - number of neighbors [int]
     - positions of neighboring atoms (relative to the given atom) [np array of np arrays containing 3 floats]
     - a np array with elements being dicts containing "distance" and "index" of neighboring atoms [np array of dictionaries]
    
    """
    # Get atom_type and aoi (atom-of-interest)
    aoi = system[ind]
    if atom_type == "*Insert Symbol Here*":
        atom_type = system[ind].symbol
    
    # Create an array with the distance between the given atom of the correct type, and all other atoms (and save the ids) using dictionaries
    distList = np.array([{"distance": atomDist(atom, aoi), "index": atom.index} 
                for atom in system if atom.symbol == atom_type and atom.index != ind])
    
    # Find the lowest distance
    lowestDist = 10**9
    for e in distList:
        if e["distance"] < lowestDist:
            lowestDist = e["distance"]
    
    # Now find the nearest neighbors, defined as every atom with a distance of less than 115% the lowest value
    for i in reversed(range(len(distList))):
        if distList[i]["distance"] > lowestDist*1.15:
            distList = np.delete(distList, i)
    
    # Create list with the relative distance (relative to aoi) of all neighboring atoms
    positions = np.array([system[e["index"]].position - aoi.position for e in distList])
    
    # Count the elements left in the distList, and return this number
    n_neighbors = len(distList)
    
    return n_neighbors, positions, distList


# Second function, cylinder (wide/tall)
def NNs_Cylinder(system, ind, shape = "wide", atom_type = "*Insert Symbol Here*"):
    """
    This function calculates the distances from the given atom, and all other atoms in the system. It then takes the shortest distance and uses it to cylinder around the given atom, with either radius or height being 115% of this distance. Everything within cylinder is counted as a neighboring atom.
    
    Use shape = "wide" for a wide cylinder, and shape = "tall" for a tall cylinder.
    
    
    Input:
     - ase atoms system [Atoms]
     - index [int]
     - shape "wide"/"tall" [string]
     - atomic symbol/type [string]
    
    Output:
     - number of neighbors [int]
     - positions of neighboring atoms (relative to the given atom) [np array of np arrays containing 3 floats]
     - a np array with elements being dicts containing "distance" and "index" of neighboring atoms [np array of dictionaries]
    
    """
    # Input sanitation
    shape = shape.lower()
    
    # Get atom_type and aoi (atom-of-interest), as well as atom radius
    atomR = 1 #PLACEHOLDER
    aoi = system[ind]
    if atom_type == "*Insert Symbol Here*":
        atom_type = system[ind].symbol
    
    # Create an array with the distance between the given atom of the correct type, and all other atoms (and save the ids) using dictionaries
    if shape == "wide":
        distList = np.array([{"distance": atomDist(atom, aoi), "index": atom.index} 
                             for atom in system 
                             if atom.symbol == atom_type  # Only count the correct type
                             and atom.index != ind  # Do not count self
                             and (aoi.position[2] - atomR < atom.position[2])  # Lower bound
                             and (aoi.position[2] + atomR > atom.position[2])]) # Upper bound
    elif shape == "tall":
        distList = np.array([{"distance": atomDist(atom, aoi), "index": atom.index} 
                             for atom in system 
                             if atom.symbol == atom_type  # Only count the correct type
                             and atom.index != ind  # Do not count self
                             and sqrt((aoi.position - atom.position)[0]**2 + (aoi.position - atom.position)[1]**2) < atomR]) # Check whether the atom below or above is more than atomR away on the x- and y-axes
    
    # Find the lowest distance
    lowestDist = 10**9
    for e in distList:
        if e["distance"] < lowestDist:
            lowestDist = e["distance"]
    
    # Now find the nearest neighbors, defined as every atom with a distance of less than 115% the lowest value
    for i in reversed(range(len(distList))):
        if distList[i]["distance"] > lowestDist*1.10:
            distList = np.delete(distList, i)
    
    # Create list with the relative distance (relative to aoi) of all neighboring atoms
    positions = np.array([system[e["index"]].position - aoi.position for e in distList])
    
    # Count the elements left in the distList, and return this number
    n_neighbors = len(distList)
    
    return n_neighbors, positions, distList


def id_array(NN_cyl_output: dict) -> array:
    """
    This function is a helper function to pack out the output of the NN_cylinder and return it in a simple list type
    it takes a NN_cyl_output and cleans the output
    """
    idArray = []
    for i in range(0,NN_cyl_output[0]):
        idArray.append(NN_cyl_output[2][i]['index'])
    return idArray
    

def get_oposite_indeces_S(system, index):
    """
    A function to get the circle of neigbohrs on the other side of the sample
    Parameters:
    -----------
     - `system` [ASE] -> the ASE system under consideration
     - `index` [int] -> the index that we consider

    Returns:
    --------
     - `cylNeigbohrsIndeces` [list] -> an array of the indeces that have all the neighbohrs on the opposite side.

    """

    # Recieving the opposite neigbbohrs index
    cylNeigbohrsTall = NNs_Cylinder(system, index, shape = "tall", atom_type = "S")
    cylNeigbohrsIndecesTall = cylNeigbohrsTall[2][0]['index']

    # Getting the neighbors around the opposite site neigbohrs indeces
    cylNeigbohrsIndeces = id_array(NNs_Cylinder(system, cylNeigbohrsIndecesTall, shape = "wide", atom_type = "S"))
    cylNeigbohrsIndeces.append(cylNeigbohrsIndecesTall)

    return cylNeigbohrsIndeces


def get_removed_indeces(system, theIndex) -> array[int]:
    """ get_removed_indeces - takes a system and an index and get all the neighborhs that exist around the atom with index and on the opposite site
    as well as all the neigbohrs to the oposite site.

    Parameters:
    -----------
     - `system` [ASE-object]: the ASE system under consideration
     - `theIndex` [int]: the index of the atom we are considering

    Returns:
    --------
     - `fullList` [array]: The list of all the integers
    
    """

    # Getting the Mo neigbohrs
    allData = NNs_Sphere(system, theIndex, "Mo")
    allDataReduced = allData[2]

    MoArray = []
    for MoAtom in range(0,len(allDataReduced)):
        MoArray.append(allData[2][MoAtom]['index'])

    # Getting the neigbohrs of the MOs
    fullList = []
    for MoAtom in MoArray:
        allData2 = NNs_Sphere(system, MoAtom, "S")
        allData2Reduced = allData2[2]
        for k in range(0,len(allData2Reduced)):
            fullList.append(allData2Reduced[k]['index'])
    
    # remove duplicates
    fullList = list(set(fullList))
    fullList.remove(theIndex)

    return fullList


def make_list_of_lists(indexList) -> array[array]:
    ''' make_list_of_lists - takes and array and make smaller combs of it,
    so for example [3,2,1] will be
    [   [1],
        [2,1],
        [3,2,1]
    ]

    Parameters:
    ----------
     - `indexList` [list]: The list we want expanded

    Returns:
    --------
     - `listOfLists`[array[array]]: see example
    '''

    listOfLists = []
    for i in range(1,len(indexList)+1):
        useList = indexList.copy()
        listOfLists.append(useList[0:i])
    return listOfLists



def get_all_combs_3(system, index, degenerate = False):
    """ get_all_combs_3 - A function which get all possibilitets of 3 of the nearest neigbohrs made with sphere_indexes function
    and uniques them such that no identical job exists

    Parameters:
    -----------
     - `system` [ASE-object]: the ASE system under consideration
     - `index` [int]: the index of the atom we are considering

    Returns:
    --------
     - `listOfLists3` [array[array]]: A list with all unqiue combinations of length 3 of the neigbhohrs removed
    """

    #Getting the neares neigbohrs and make all possible unqiue combinations
    ids = sphere_indexes(system, index)
    listOfLists3 = list()
    for i in ids:
        for j in ids:
            for k in ids:
                if(i!=j and i!=k and j!=k):
                    listen = [i,j,k]
                    listen.sort()
                    listOfLists3.append(listen)
    if(degenerate == False):
        listOfLists3 = make_unique(listOfLists3)
    return listOfLists3

def get_all_combs_2(system, index, degenerate = False):
    """ get_all_combs_2 - A function which get all possibilitets of 2 of the nearest neigbohrs made with sphere_indexes function
    and uniques them such that no identical job exists

    Parameters:
    -----------
     - `system` [ASE-object]: the ASE system under consideration
     - `index` [int]: the index of the atom we are considering

    Returns:
    --------
    - `listOfLists2` [array[array]]: A list with all unqiue combinations of length 2 of the neigbhohrs removed
    """

    #Getting the neares neigbohrs and make all possible unqiue combinations
    ids = sphere_indexes(system, index)
    listOfLists2 = list()
    for i in ids:
        for j in ids:
            if(i!=j):
                listen = [i,j]
                listen.sort()
                listOfLists2.append(listen)
    if(degenerate == False):
        listOfLists2 = make_unique(listOfLists2)
    return listOfLists2
                
def get_all_combs_1(system, index, degenerate = False): 
    """ get_all_combs_1 - A function which get all possibilitets of 1 of the nearest neigbohrs made with sphere_indexes function
    and uniques them such that no identical job exists

    Parameters:
    -----------
     - `system` [ASE-object]: the ASE system under consideration
     - `index` [int]: the index of the atom we are considering

    Returns:
    --------
     - `listOfLists1` [array[array]]: A list with all unqiue combinations of length 1 of the neigbhohrs removed
    """
    
    #Getting the neares neigbohrs and make all possible unqiue combinations
    ids = sphere_indexes(system, index)
    listOfLists1 = list()
    for i in ids:
        listen = [i]
        listen.sort()
        listOfLists1.append(listen)
    if(degenerate == False):
        listOfLists1 = make_unique(listOfLists1)
    return listOfLists1

def make_unique(listOfLists):
    """ make_unique - A function which takes an array of arrays and removes make it so only unqiue occurences exist
    
    Example:
    --------
    [[1,2,3],[1,3,2],[1,2,3]]   --f-->    [[1,2,3],[1,3,2]] 
    """
    # Unique the combinations
    uniqueList = []
    for j in listOfLists:
        if(j not in uniqueList):
            uniqueList.append(j)
    return uniqueList

def get_all_combs_tot(system, index, degenerate = False):
    """ get_all_combs_tot - takes a system and the index of the atom we consider and makes a List containing a list of
    all the atoms we can remove of the combinations of 1, 2, 3 removals

    Parameters:
    -----------
     - `system` [ASE-object]: the ASE system under consideration
     - `index` [int]: the index of the atom we are considering

    Returns:
    --------
     - `listTot` [array[array]]: A list with all unqiue combinations of length 1, 2 and 3 neighbors to be removed
    """

    # collecting lists and combine them
    list1 = get_all_combs_1(system, index, degenerate)
    list2 = get_all_combs_2(system, index, degenerate)
    list3 = get_all_combs_3(system, index, degenerate)
    listTot = list1 + list2 + list3
    return listTot
    
def sphere_indexes(system, index) -> list:
    """ sphere_indexes - A quick helper function to get out the nearest neighbohrs of the spheres indexes
    it is like calling a NNs_Sphere but more simple
    """

    sphereOut = NNs_Sphere(system, index, "S")
    sphereOut = sphereOut[2]
    sphereOutArray = list()
    for element in sphereOut:
        sphereOutArray.append(element['index'])
    
    return sphereOutArray




