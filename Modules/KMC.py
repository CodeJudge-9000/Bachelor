
# Initialization
import FingerPrints
from atomic_annihilator import *
import numpy as np
import math as m
from math import sqrt
from random import uniform
from random import randint
from copy import copy


class KMC:
    def __init__(self, tempSize, TDlib, fingerPrint, kinetic_E, electron_dose):
        self.create_system(tempSize)
        self.dose = electron_dose # number of electrons/(Å**2 * s)
        self.total_sim_time = 0 # Some time unit. Figure this one out later
        self.current_sim_time = 0 # Same, figure out the time later
        self.rate_constant = self.get_rate_constant() # Calculate the rate constant for the system
        self.energy_cutoff_S = self.get_energy_cutoff() # Calculate the energy cutoff for our atom (S)
        self.b_cutoff_S = self.get_b_cutoff() # Calculate the b-cutoff for our atom (S)
        self.TDlib = TDlib # TD energies in eV
        self.fingerPrint = fingerPrint # Name of the fingerPrint used
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV)
        self.missingTDs = [] # Add the missing fingerprints to this list
    


    def set_electron_energy(self, kinetic_E):
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV)
        return
    
    def set_electron_dose(self, electron_dose):
        self.dose = electron_dose
        return


    def run(self):
        self.current_sim_time = 0
        
        return
    


    def simulate_electron(self):
        # Choose which electron to interact with
        sideLen = len(self.grid_S[-1])

        a1, a2 = randint(0, sideLen-1), randint(0, sideLen-1)

        # Figure out the interaction distance (b-value) of the electron and atom
        b_cutoff = self.get_b_cutoff()


        p1 = uniform(0, b_cutoff**2)
        b = m.sqrt(p1)

        # Something something energy transferred depending on distance


        return
    


    def create_system(self, tempSize):
        # Create S grid
        self.grid_S = np.ones((3, tempSize, tempSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)
        self.grid_S[1][:][:] = False # Set the middle-layer to be false

        # Create Mo grid
        self.grid_Mo = np.ones((tempSize-1, tempSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

        return
    


    def get_transferred_energy(self):

        return



    def get_rate_constant(self):
        # Determine how many atoms can get hit (In this case the amount of atoms in the bottom layer)
        numberS = np.sum(self.grid_S[-1])

        # Determine the rate constant
        rate_constant = self.dose * m.pi * self.get_b_cutoff()**2 * numberS # 1/s

        return rate_constant
        
    
    def get_b_cutoff(self):
        # Use the highest TD value to find the b cutoff

        return

    def get_energy_cutoff(self):

        return


    def time_step(self):
        # Calculate how long passed
        timePassed = -1 * m.log(uniform(0.0000000001, 1.0))/self.rate_constant # Change in seconds

        # Total time ran update
        self.total_sim_time += timePassed

        # This simulation cycle update
        self.current_sim_time += timePassed

        return



    def get_TD(self, index):
        """
        Using an id, calculate the fingerprint of the atom and return its TD value, using the TD library. If no TD value exists for the given fingerprint, log this (add to the missingTDs list) and return False.
        """
        # First get the fingerprint for the corresponding atom id
        finger_method = getattr(FingerPrints, self.fingerPrint)
        finger = finger_method(self.system, index)
        
        ### BUGTESTING ###
        if finger[0] == 4 or finger[0] == 6:
            print(index, finger)
        ### BUGTESTING ###
        
        # Now run through the pandas dataframe, and check if there are any corresponding value
        if len(self.TDlib[self.TDlib["finger"] == str(finger)]) == 0:
            # If there are none, add the fingerPrint to missingTDs (if it is not there already) and return True
            if finger not in self.missingTDs:
                self.missingTDs.append(finger)
            return False
        
        # If there are at least one corresponding TD value, take the average of all the values and return the value
        else:
            return self.TDlib[self.TDlib["finger"] == str(finger)].mean()["Td"]
    
    

    def get_missing_TDs(self):
        """
        Returns a list of all the fingerprints missing a TD value.
        """
        return self.missingTDs
    


    def clear_missing_TDs(self):
        "Clears the missing TD list."
        self.missingTDs = []
        return






# 2D distance function
def Dist_2D(pos1, pos2):
    d1 = pos1[0] - pos2[0] # The x-axis
    d2 = pos1[1] - pos2[1] # The y-axis
    return sqrt(d1**2 + d2**2)


class KMC_OLD:
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
        
        tempmissingTDs = copy(self.missingTDs)
        for i in range(n):
            self.simulate_electron()
        
        return tempmissingTDs == self.missingTDs
        
    
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
            print("Removed one")
            remove_atoms(self.system, atomRemoveIndex = lowestIndex, relax = False, overwriteCalc = True)
        
        return True
    
    def cross_section(self):
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
        
    def cross_dummy(self, symbol):
        """
        Dummy-function to imitate calculating the interaction cross-section of an atom.
        """
        if symbol == "S":
            return 2
        else:
            return 5
    
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
        x = uniform(0, xLen) # TODO
        y = uniform(0, yLen) # TODO
        point = (x, y)
        
        # Get the center of mass' third coordinates
        comZ = self.system.get_center_of_mass()[2]
        
        # Check for all atoms (that are S) below the center of mass if any cross-sections are intersected
        interSect = [(atom.index, Dist_2D(atom.position, point))
                       for atom in self.system if atom.position[2] < comZ
                       and atom.symbol == "S"
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
        finger_method = getattr(FingerPrints, self.fingerPrint)
        finger = finger_method(self.system, index)
        
        ### BUGTESTING ###
        if finger[0] == 4 or finger[0] == 6:
            print(index, finger)
        ### BUGTESTING ###
        
        # Now run through the pandas dataframe, and check if there are any corresponding value
        if len(self.TDlib[self.TDlib["finger"] == str(finger)]) == 0:
            # If there are none, add the fingerPrint to missingTDs (if it is not there already) and return True
            if finger not in self.missingTDs:
                self.missingTDs.append(finger)
            return False
        
        # If there are at least one corresponding TD value, take the average of all the values and return the value
        else:
            return self.TDlib[self.TDlib["finger"] == str(finger)].mean()["Td"]
    
    
    def get_missing_TDs(self):
        """
        Returns a list of all the fingerprints missing a TD value.
        """
        return self.missingTDs
    
    def clear_missing_TDs(self):
        "Clears the missing TD list."
        self.missingTDs = []
        return
