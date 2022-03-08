import ase.build
import numpy as np


def RemovBelowLine(fCoord, sCoord, system):
    """
    Remove particles under a line drawn between two particles (ignores depth)
    fCoord should be the coordinate of the first particle
    sCoord should be the coordinate of the second particle (to the right of the first particle)
    
    Both should be numpy arrays"""
    
    
    # Get the x- and y-coordinates
    xC = np.array([fCoord[0], sCoord[0]])
    yC = np.array([fCoord[1], sCoord[1]])
    
    # Fit to a linear function, and add a bit to the 0'th degree parameter
    pParams = np.polyfit(xC, yC, 1) + np.array([0, 0.0001])

    # Remove atoms using list comprehension
    del system[[atom.index for atom in system if (atom.position[1]<=atom.position[0]*pParams[0]+pParams[1])]]

    return    

def lBelowRemover(fIndx, sIndx, system):
    """
    Remove particles under a line drawn between two particles for the given system (ignores depth)
    fIndx should be the index of the first particle (the particle to draw from)
    sIndx should be the index of the second particle (the particle to draw to)"""
    
    # Get the two particle positions
    fCoord = system.get_positions()[fIndx]
    sCoord = system.get_positions()[sIndx]
    
    # Get the x- and y-coordinates
    xC = np.array([fCoord[0], sCoord[0]])
    yC = np.array([fCoord[1], sCoord[1]])
    
    # Fit to a linear function, and add a bit to the 0'th degree parameter
    pParams = np.polyfit(xC, yC, 1) + np.array([0, 0.0001])

    del system[[atom.index for atom in system if (atom.position[1]<=atom.position[0]*pParams[0]+pParams[1])]]

    return

def vBelowRemover(i1, i2, i3, system):
    """Removes everything below a line drawn from particle with index i1, to i2, to i3"""

    # Get the coordinates for each atom
    Coord1 = system.get_positions()[i1]
    Coord2 = system.get_positions()[i2]
    Coord3 = system.get_positions()[i3]
    
    # Remove under lines using RemovBelowLine twice
    RemovBelowLine(Coord1, Coord2, system)
    RemovBelowLine(Coord2, Coord3, system)

    return