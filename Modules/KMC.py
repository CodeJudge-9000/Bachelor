
# Initialization
from FingerPrints import *
import numpy as np
from math import sqrt
from random import randrange

# 2D distance function
def 2D_Dist(pos1, pos2):
    d1 = pos1[0] - pos2[0] # The x-axis
    d2 = pos1[1] - pos2[1] # The y-axis
    return sqrt(d1**2 + d2**2)


class KMC:
    def __init__(self, system, TDlib, fingerPrint):
        self.system = system
        self.TDlib = TDlib
        self.fingerPrint = fingerPrint # Name of the fingerPrint used
        self.missingTDs = [] # Add the missing fingerprints to this list
    
    def run(n):
        """
        Simulates n electrons, and return True or False depending on whether any TD values were missing. If no TD value was missing, return False, otherwise return True.
        """
        tF = 0
        for i in range(n):
            tF += self.simulate_electron() # simulate_electron() returns False (0) of no TD value was missing, and True (1) if a value was missing
        
        return bool(tF)
        
    
    def simulate_electron():
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
        
        #
        # 1. Determine where to hit structure
        # 2. Do something with the cross-section of atoms nearby (perhaps a search?), to determine which atom(s) are affected
        # 2.1 Check whether a TD value exists for this/these atom(s). If not, return 1
        # 3. Determine how much energy is transferred, if energy is transferred
        # - Use mr. PHD's distribution for this
        # 4. If above threshold value, remove the atom
        # 4.1 Maybe implement some other behaviour, such as removing neighboring atoms
        # 4.2 Also perhaps add the energy to the system in case it does not exceed the TD
        #
        # Return 0 at the end (if the TD value for the atom exists)
        #
    
    def cross_section():
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
        
    def cross_dummy():
        """
        Dummy-function to imitate calculating the interaction cross-section of an atom.
        """
        if system[index].symbol = "S":
            return 1
        else:
            return 1.5
    
    def determine_hit():
        """
        Determines where an electron enters the structures, and determines whether it'll intersect any cross-section for an atom. It then returns the ID of this atom if True, otherwise it returns False.
        """
        # Get cell dimensions
        cell = self.system.get_cell() # get_cell returns the lengths of the cell, as 3 vectors stretching out from the point (0, 0, 0)
        xLen = cell[0][0]
        yLen = cell[1][1]
        
        # Choose a random point (x, y) within the cell
        x = randrange(0, xLen) # TODO
        y = randrange(0, yLen) # TODO
        
        # Get the center of mass' third coordinates
        comZ = system.get_center_of_mass()[2]
        
        # Check for all atoms below the center of mass if any cross-sections are intersected
        [atom.index for atom in system if atom.position[2] < comZ and 2D_Dist(atom.position, )] # TODO
        
        
    def get_TD(index):
        """
        Using an id, calculate the fingerprint of the atom and return its TD value, using the TD library. If no TD value exists for the given fingerprint, log this (add to the missingTDs list) and return False.
        """
        # First get the fingerprint for the corresponding atom id
        finger_method = getattr(FingerPrints, fingerPrint)
        finger = finger_method(index)
        
        # Now run through the pandas dataframe, and check if there are any corresponding value
        if len(self.TDlib[self.TDlib["finger"] == finger]) == 0:
            # If there are none, add the fingerPrint to missingTDs and return True
            self.missingTDs.append(finger)
            return False
        
        # If there are at least one corresponding TD value, take the average of all the values and return the value
        else:
            return self.TDlib[self.TDlib["finger"] == finger].mean()["Td"]
    
    
    def get_missing_TDs():
        """
        Returns a list of all the fingerprints missing a TD value.
        """
        return self.missingTDs
    
    
