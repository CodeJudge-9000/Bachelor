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
        if distList[i]["distance"] > lowestDist*1.15:
            distList = np.delete(distList, i)
    
    # Create list with the relative distance (relative to aoi) of all neighboring atoms
    positions = np.array([system[e["index"]].position - aoi.position for e in distList])
    
    # Count the elements left in the distList, and return this number
    n_neighbors = len(distList)
    
    return n_neighbors, positions, distList

def id_array(NN_cyl_output):
    idArray = []
    for i in range(0,NN_cyl_output[0]):
        idArray.append(NN_cyl_output[2][i]['index'])
    return idArray
    
def get_oposite_indeces_S(sys, ind):
    cylNeigbohrsTall = NNs_Cylinder(system = sys, ind = ind, shape = "tall", atom_type = "S")
    cylNeigbohrsIndecesTall = cylNeigbohrsTall[2][0]['index']
    cylNeigbohrs = NNs_Cylinder(system = sys, ind = cylNeigbohrsIndecesTall, shape = "wide", atom_type = "S")
    cylNeigbohrsIndeces = id_array(cylNeigbohrs)
    cylNeigbohrsIndeces.append(cylNeigbohrsIndecesTall)
    return(cylNeigbohrsIndeces)

def get_removed_indeces(sys, theIndex):

    allData = NNs_Sphere(system = sys, ind = theIndex, atom_type = "Mo")
    allDataReduced = allData[2]

    MoArray = []
    for i in range(0,len(allDataReduced)):
        MoArray.append(allData[2][i]['index'])

    fullList = []
    for j in MoArray:
        allData2 = NNs_Sphere(system = sys, ind = j, atom_type = "S")
        allData2Reduced = allData2[2]
        for k in range(0,len(allData2Reduced)):
            fullList.append(allData2Reduced[k]['index'])
    fullList = list(set(fullList))
    fullList.remove(theIndex)

    return fullList


def make_list_of_lists(indexList):
    listOfLists = []
    for i in range(1,len(indexList)):
        useList = indexList.copy()
        listOfLists.append(useList[0:i])
    return listOfLists


