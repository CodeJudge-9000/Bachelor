
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
        # First create the system
        self.create_system(tempSize)

        # Define a bunch of internal parameters
        self.dose = electron_dose # number of electrons/(Å**2 * s)
        self.total_sim_time = 0 # Some time unit. Figure this one out later
        self.TDlib = TDlib # TD energies in eV
        self.fingerPrint = fingerPrint # Name of the fingerPrint used
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV) <- Units are *important*

        # Define some constants
        self.m_e = 9.1093837015*10**(-31) # Electron mass in kg
        self.m_S = 32.06 * 1.660540200*10**(-27) # Sulfur mass in kg
        self.m_Mo = 95.95 * 1.660540200*10**(-27) # Molybdenum mass in kg
        self.speed_of_light_si = 2.998*10**8 # 'c' in m/s
        self.coulomb_k_si = 8.9875517923*10**9 # The coulomb constant in SI units

        # Calculate some of the constant values of the system
        self.rate_constant_S = self.get_rate_constant("S") # Calculate the rate constant for the system
        self.relativistic_electron_mass = self.get_relativistic_electron_mass() # The relativistic electron mass in kg
        self.electron_velocity = self.get_electron_velocity() # Electron velocity in m/s
        self.energy_cutoff_S = self.get_energy_cutoff("S") # Calculate the energy cutoff for our atom type (S)
        self.b_cutoff_S = self.get_b_cutoff("S") # Calculate the b-cutoff for our atom type (S)
        self.a_S = self.a("S") # Calculate 'a'

        # Create the initial (and empty) missing fingerprints list
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
    


    def simulate_electron(self): # Fix this thing
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
    


    def get_transferred_energy(self, b, atomSymb):
        """Calculates and returns the transferred energy in eV"""
        # First calculate the momentum
        p_trans = (2*self.relativistic_electron_mass * self.electron_velocity) / m.sqrt(b**2 * self.a_S**2 + 1)

        # Then find other required parameters
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo

        # Then use the momentum to calculate the energy
        E_T = 6.241509125*10^18*p_trans**2/(2*m_n)

        return E_T



    def get_rate_constant(self, atomSymb):
        """Calculates and returns the rate constant of the whole system for ONE type of atom in 1/s"""
        # Determine how many atoms can get hit (In this case the amount of atoms in the bottom layer)
        numberS = np.sum(self.grid_S[-1])

        # Determine the rate constant
        rate_constant = self.dose * m.pi * self.get_b_cutoff(atomSymb)**2 * numberS # 1/s

        return rate_constant

    def get_relativistic_electron_velocity(self):
        """Calculates and returns the speed of the electrons, considering relativistic effects"""
        v_rela = 2.998*10**8*m.sqrt(1 - 1/(self.electronKin/(5.1098895*10**5) + 1)**2)
        return v_rela
    
    def get_electron_velocity(self):
        """Calculates and returns the classical electron velocity in m/s"""
        v_rest = m.sqrt(2*self.electronKin/(self.m_e*6.241509125*10**(-18)))

        return v_rest
    
    def get_relativistic_electron_mass(self):
        """Calculates the relativistic mass of the electron and returns it in kg"""
        v = self.get_electron_velocity()
        m_r = (self.m_e) / (1 - (v/self.speed_of_light_si)**2)
        return m_r

    def a(self, atomSymb):
        """Calculates the 'a' parameter, and returns it in 1/m"""
        # First get the variables in order
        v_0 = self.get_electron_velocity()

        if atomSymb == "S":
            m_n = self.m_S
            Q = 16
        elif atomSymb == "Mo":
            m_n = self.m_Mo
            Q = 42

        # Now calculate a
        a = (v_0**2 * self.m_e) / (self.coulomb_k_si * (-1) * Q * (self.m_e/m_n + 1)**3)

        return a
    
    def get_b_cutoff(self, atomSymb):
        """Calculates and returns the cutoff value for b in Å (angstrom)"""
        # Find the lowest TD value, as to find the b cutoff (as E_T ~ 1/b**2)
        E_min = self.TDlib["Td"].min() * 1.05

        # Now get the mass of the atomic nucleus of the corresponding atom
        # The following should (for maximum compatibility) by some library but for now it's just some if-else statements
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo

        # Get the relativistic mass of our electrons, as well as their velocity
        m_r = self.get_relativistic_electron_mass()
        v_0 = self.get_electron_velocity()

        # Get 'a'
        a = self.a(atomSymb)

        # Calculate the cutoff value for b
        b_cutoff = 10**(-10)/a * m.sqrt((2*m_r*v_0)**2 / (6.241509125*10**(-18)*E_min*2*m_n))

        # Convert it to Å and return it
        b_cutoff = b_cutoff * 10**(-10)

        return b_cutoff

    def get_energy_cutoff(self, atomSymb):
        """Calculates and returns the energy cutoff (maximum) in eV"""
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo

        E_max = (2* self.electronKin * (self.electronKin + 2*5.1098895*10**5)) / (6.241509125*10**18 * m_n * self.speed_of_light_si**2)

        return E_max


    def time_step(self):
        """Increases the total sim time as well as the current sim time using the total rate constant for the system"""
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
        """Clears the missing TD list."""
        self.missingTDs = []
        return

