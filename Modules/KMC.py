
# Initialization
from FingerPrints import *
from atomic_annihilator import *
import numpy as np
from math import sqrt
from random import randrange

# 2D distance function
def Dist_2D(pos1, pos2):
    d1 = pos1[0] - pos2[0] # The x-axis
    d2 = pos1[1] - pos2[1] # The y-axis
    return sqrt(d1**2 + d2**2)


class KMC:
    def __init__(self, system, TDlib, fingerPrint, eE):
        self.system = system.copy() # Create a copy of the ase system
        self.TDlib = TDlib # TD energies in eV
        self.fingerPrint = fingerPrint # Name of the fingerPrint used
        self.electronE = eE # Same unit as TD values (eV?)
        self.missingTDs = [] # Add the missing fingerprints to this list
    
    def set_electron_energy(self, eE):
        self.electronE = eE # eV
        return
    
    def run(self, n):
        """
        Simulates n electrons, and return True or False depending on whether any TD values were missing. If no TD value was missing, return True, otherwise return False.
        """
        tF = True
        for i in range(n):
            simBool = self.simulate_electron()
            if simBool == False:
                tF = False
        
        return tF
        
    
    def simulate_electron(self):
        
        #
        # Note: X - Done, ? - Dummy/testing function implemented, M - Not implemented
        #
        # 1.X Determine where to hit structure (Part of determine_hit. It chooses a x and y coordinate seperately from an uniform distribution)
        # 2.X Do something with the cross-section of atoms nearby (perhaps a search?), to determine which atom(s) are affected
        # 3.X Check whether a TD value exists for this/these atom(s). If not, return False
        # 4.? Determine how much energy is transferred, if energy is transferred
        # - Use mr. PHD's distribution for this
        # 5.X If above threshold value, remove the atom
        # 5.1.M Maybe implement some other behaviour, such as removing neighboring atoms
        # 5.2.M Also perhaps add the energy to the system in case it does not exceed the TD
        #
        # X Return True at the end (if the TD value for the atom exists), otherwise Return False at some point
        #
        
        # Determine whether anything is hit by the electron
        interSect = self.determine_hit()
        
        if interSect == False:
            return
        
        # Now that we have one (or more) hits, we choose the closest intersect
        lowestDist = 99999
        lowestIndex = None
        for e in interSect:
            if e[1] < lowestDist:
                lowestDist = e[1]
                lowestIndex = e[0]
        
        # Get the TD value for the closest intersect
        TD = self.get_TD(lowestIndex)
        if TD == False:
            return False
        
        # Now we figure out how much energy is transferred
        energyTransfer = self.dummy_energy_transferred()
        
        # If this exceeds the TD value for the atom, remove it (use the atomRemover function as to remove a list easily, and possibly keep track of indices in the future)
        if  TD <= energyTransfer:
            remove_atoms(self.system, atomRemoveIndex = lowestIndex, relax = False, overwriteCalc = True)
        
        return
    
    def cross_section(self):
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
        
    def cross_dummy(self, index):
        """
        Dummy-function to imitate calculating the interaction cross-section of an atom.
        """
        if self.system[index].symbol = "S":
            return 1
        else:
            return 1.5
    
    def energy_transferred(self):
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
    
    def dummy_energy_transferred(self):
        """
        Dummy-function to imitate calculating the energy transferred from an electron to an atom.
        """
        return self.electronE
        
    
    def determine_hit(self):
        """
        Determines where an electron enters the structures, and determines whether it'll intersect any cross-section for any atoms. It then returns a list containing the IDs and distances of these atom if True (tuples in a list), otherwise it returns False.
        
        The returned list is in the following format:
            [(index, distance_to_atom), (index, distance_to_atom), ...]
        """
        # Get cell dimensions
        cell = self.system.get_cell() # get_cell returns the lengths of the cell, as 3 vectors stretching out from the point (0, 0, 0)
        xLen = cell[0][0]
        yLen = cell[1][1]
        
        # Choose a random point (x, y) within the cell
        x = randrange(0, xLen) # TODO
        y = randrange(0, yLen) # TODO
        point = (x, y)
        
        # Get the center of mass' third coordinates
        comZ = system.get_center_of_mass()[2]
        
        # Check for all atoms below the center of mass if any cross-sections are intersected
        interSect = [(atom.index, Dist_2D(atom.position, point))
                       for atom in system if atom.position[2] < comZ
                       and Dist_2D(atom.position, point) <= self.cross_dummy(atom.symbol)] # REPLACE cross_dummy() WITH cross_section() WHEN POSSIBLE
        # If no intercepts were found, return false
        if len(interSect) == 0:
            return False
        
        # Else, return a list of tuples, with each tuple containing the atomic index, and the distance from the electron to the atom (usually just return a tuple, but for futureproofing we return a list)
        else:
            return interSect
        
    def get_TD(self, index):
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
    
    
    def get_missing_TDs(self):
        """
        Returns a list of all the fingerprints missing a TD value.
        """
        return self.missingTDs
    
    
