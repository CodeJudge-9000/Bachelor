
# Initialization
from ase import Atoms
from ase.io.trajectory import TrajectoryWriter
import warnings
import numpy as np
import math as m
from random import uniform
from random import randint


class KMC:
    def __init__(self, gridSize, structureType, TDlib, fingerPrint, kinetic_E, electron_dose_rate, overwriteTD = False):
        # First create the system
        self.create_system(gridSize, structureType)

        # Define a bunch of internal parameters
        self.overwriteTD = overwriteTD
        self.squareSize = gridSize
        self.structType = structureType
        self.dose_rate = electron_dose_rate # number of electrons/(Å**2 * s)
        self.total_sim_time = 0 # Some time unit. Figure this one out later
        self.fingerPrint = fingerPrint
        self.TDlib = TDlib # TD energies in eV
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV) <- Units are *important*
        self.gridStack = [(self.grid_S.copy(), self.grid_Mo.copy(), self.grid_Removed.copy(), self.total_sim_time, 0)] # list used to store the two grids, as well as the 
        self.S_init = np.sum(self.grid_S)

        # Define some constants
        self.m_e = 9.1093837015*10**(-31) # Electron mass in kg
        self.m_S = 32.06 * 1.660540200*10**(-27) # Sulfur mass in kg
        self.m_Mo = 95.95 * 1.660540200*10**(-27) # Molybdenum mass in kg
        self.speed_of_light_si = 2.998*10**8 # 'c' in m/s
        self.coulomb_k_si = 2.307077515*10**(-28) # The coulomb constant times e^2, given in SI units (kg*m^3)/(s^2)
        self.k_B_si = 1.380649*10**(-23) # The boltzmann constant given in J/K

        # Calculate some of the constant values of the system
        self.b_cutoff_mean_S = self.get_b_cutoff("S",0) # Calculate the b-cutoff for our atom type (S)
        self.rate_constant_S = self.get_rate_constant("S") # Calculate the rate constant for the system
        self.relativistic_electron_mass = self.get_relativistic_electron_mass() # The relativistic electron mass in kg
        self.electron_velocity = self.get_relativistic_electron_velocity() # Electron velocity in m/s
        self.energy_cutoff_mean_S = self.get_energy_cutoff("S",0) # Calculate the energy cutoff for our atom type (S)
        self.a_S = self.a("S") # Calculate 'a'
        self.m_r_eS = self.get_reduced_mass("S")

        # Create the initial (and empty) missing fingerprints list
        self.missingTDs = [] # Add the missing fingerprints to this list
    


    def set_electron_energy(self, kinetic_E):
        """Updates the electron energy"""
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV)

        # Remember to update the dependent functions
        self.electron_velocity = self.get_electron_velocity()
        self.energy_cutoff_mean_S = self.get_energy_cutoff("S",0)
        self.relativistic_electron_mass = self.get_relativistic_electron_mass()
        self.a_S = self.a("S")
        self.b_cutoff_mean_S = self.get_b_cutoff("S",0)
        self.rate_constant_S = self.get_rate_constant("S")
        return
    
    def set_electron_dose_rate(self, electron_dose_rate):
        """Updates the electron dose rate"""
        self.dose_rate = electron_dose_rate

        # Remember to update the dependent functions
        self.rate_constant_S = self.get_rate_constant("S")
        return

    def run(self, iterN, feedBack = False, sample = True, runToLimit = False, useT = True):
        self.current_sim_time = 0
        iterations = 0
        
        missedElectrons = 0
        
        if runToLimit == True:
            sumAtoms = self.S_init
            while sumAtoms > self.S_init*0.90:
                outcome = self.simulate_electron(sample, useT)
                if  outcome == False:
                    missedElectrons += 1
                else:
                    self.update_rate_constant_S()
                    self.time_step()
                    if outcome == True:
                        sumAtoms -= 1
                    else:
                        missedElectrons += 1
                    
                self.gridStack.append((self.grid_S.copy(), self.grid_Mo.copy(), self.grid_Removed.copy(), self.total_sim_time))
                iterations += 1
                
                # Condition to ensure we do not end up in an infinite loop
                termMult = 0.8
                if iterations >= self.S_init * termMult:
                    print(f"Exited loop due to termination condition iterations >= {self.S_init * termMult}")
                    break

        else:
            while iterations < iterN:
                outcome = self.simulate_electron(sample)
                if outcome == False or outcome == "interaction":
                    missedElectrons += 1
                else:
                    self.update_rate_constant_S()
                    self.time_step()
                self.gridStack.append((self.grid_S.copy(), self.grid_Mo.copy(), self.grid_Removed.copy(), self.total_sim_time, iterations))
                iterations += 1
        
        if feedBack == True:
            print(f"Number of electrons simulated: {iterations}")
            print(f"Number of missed electrons: {missedElectrons}")
            print(f"Number of atoms knocked out of structure: {iterations - missedElectrons}")
        
        return
    
    def simulate_electron(self, sample = True, useT = True): # Fixed this thing
        """Returns True if an atom was removed, and False otherwise. If an atom was interacted with, but did not get removed, this function returns the string "interaction". """
        # Choose a point in the S grid independently of the layer
        sideLen = len(self.grid_S[-1])
        a1, a2 = randint(0, sideLen-1), randint(0, sideLen-1)

        # Determine of there is an atom in any of the layers, starting from the top layer
        layer = 0
        fingerPrint = None
        while (layer < 3):
            fingerPrint = self.get_fingerPrint(layer, a1, a2)
            
            # The interaction for the top layer
            if fingerPrint != None and layer == 0:
            
                # In the case there is no atom on the other side of the structure
                if self.grid_S[2][a1,a2] == False and self.higherThanTD(fingerPrint, situation = "PT", sample = sample, useT = useT) == True:
                    # Remove atom in top-layer & add it to the bottom layer
                    self.grid_S[layer][a1,a2] = False
                    self.grid_S[2][a1,a2] = True
            
                    # Update the grid of removed atoms
                    self.grid_Removed[layer][a1,a2] = True
                    self.grid_Removed[2][a1,a2] = False
            
                    # Then update the total rate constant for the system
                    self.rate_constant_S = self.get_rate_constant("S")
                    
            
                    return True
                
                # Then in the case there is (or the transferred energy was too low) then pass (assume the electron does NOT interact with the first atom and just flies through unhindered, though in reality it may be turned around)
                else:
                    pass
            
            # Interaction for the middle-layer and bottom layer
            if fingerPrint != None and layer != 0:
                if self.higherThanTD(fingerPrint, situation = "PA", sample = sample, useT = useT) == True:
                    self.removeAtom(layer, a1, a2)
                    return True
                else:
                    pass
            # Update layer variable
            layer += 1
        
        # Check whether there was an atom in any layer, but that we missed it
        for l in [0,1,2]:
            if self.grid_S[l][a1,a2] == True:
                return "interaction"

        return False
    
    def higherThanTD(self, fingerPrint, situation, sample, useT):
        # Figure out the interaction distance (b-value) of the electron and atom
        b = m.sqrt(uniform(0, self.b_cutoff_mean_S**2)/m.pi)

        # Get a velocity for the atom
        if useT == True:
            velocity = self.get_velocity_20C("S")
        else:
            velocity = 0
        
        # Figure out how much energy is transferred & the cutoff
        E_T = self.get_transferred_energy(b, "S", velocity)
        E_cutoff = self.get_energy_cutoff("S", velocity)


        if E_T > E_cutoff:
            E_T = E_cutoff
        
        # Now check whether the transferred energy is higher than the TD value for this atom
        TD = self.get_TD(fingerPrint, situation, sample)
        
        if TD == None or E_T < TD:
            # Return False if there is no corresponding TD value or the transferred energy is lower than the TD value
            return False
        else:
            # Return True if there is a corresponding TD value, and if the transferred energy if high enough
            return True

    def removeAtom(self, layer, a1, a2):
        # Remove atom in layer
        self.grid_S[layer][a1,a2] = False

        # Also update the grid of removed atoms
        self.grid_Removed[layer][a1,a2] = True
        
        # Check if it was a corner atom, and if so remove the atom in the other layer
        if a1 + a2 == 0 or (a1 == 0 and a2 == self.squareSize-1) or (a2 == 0 and a1 == self.squareSize-1) or a1 + a2 == (self.squareSize-1)*2:
            # Remember that this only happens in the bottom layer
            if layer == 2:
                self.grid_S[0][a1,a2] = False
                self.grid_Removed[0][a1,a2] = True

        # Then update the total rate constant for the system
        self.rate_constant_S = self.get_rate_constant("S")
        return

    def reset(self):
        """Creates and updates the current grids with undamaged ones using the size of the current grids. Also resets the total simulation time and gridStack"""
        self.create_system(self.squareSize, self.structType)
        self.total_sim_time = 0
        self.gridStack = [(self.grid_S, self.grid_Mo, self.grid_Removed, self.total_sim_time, 0)]

        return

    def current_grid_to_atoms(self):
        """Converts the current grids into an ase atoms object, and returns it"""
        return self.grid_to_atoms(-1)

    def grid_to_atoms(self, stackLayer, fancy = True):
        """Converts the given gridStack grids into an ase atoms object, and returns it"""
        # Get the S and Mo grid
        grid_S = self.gridStack[stackLayer][0]
        grid_Mo = self.gridStack[stackLayer][1]
        grid_Removed = self.gridStack[stackLayer][2]

        atomS_list = []
        atomMo_list = []
        atomRemoved_list = []
        intera1 = 3.18
        intera2 = (1.59,2.754)
        interL = 1.595

        # First construct the list of S atoms
        # Go over each layer
        for L in range(grid_S.shape[0]):
            # Go over each first coordinate (aka row 'r' or a1)
            for a1 in range(grid_S.shape[1]):
                # Go over each second coordinate (aka column 'c' or a2)
                for a2 in range(grid_S.shape[2]):
                    # Begin at the first column in the first layer
                    if grid_S[L][a1,a2] == True:
                        atomS_list.append([a1*intera1 + a2*intera2[0],a2*intera2[1],(L-1)*interL])
                    if grid_Removed[L][a1,a2] == True and fancy == True:
                        atomRemoved_list.append([a1*intera1 + a2*intera2[0],a2*intera2[1],(L-1)*interL])

        # Now construct the list of Mo atoms, using much the same method
        # Go over each first coordinate (aka row 'r' or a1)
        for a1 in range(grid_Mo.shape[0]):
            # Go over each second coordinate (aka column 'c' or a2)
            for a2 in range(grid_Mo.shape[1]):
                # Begin at the first column in the first layer
                if grid_Mo[a1,a2] == True:
                    atomMo_list.append([a1*intera1 + a2*intera2[0] + 3.18,a2*intera2[1] + 1.836, 0])

        # Now that we've constructed the coordinate lists, we construct the name of the system
        sysString = f"S{len(atomS_list)}Mo{len(atomMo_list)}O{len(atomRemoved_list)}"

        # And now we finally construct the Atoms object
        system = Atoms(sysString,
                    positions=atomS_list + atomMo_list + atomRemoved_list,
                    cell=[1,1,1],
                    pbc=[0, 0, 0])

        system.center(vacuum = 1)
        
        return system
    
    def gridStack_TrajectoryWriter(self, filename):
        """Converts the current gridStack to atoms objects and writes them to a .traj file"""
        writer = TrajectoryWriter(f"{filename}.traj", "w")

        for i in range(len(self.gridStack)):
            sys = self.grid_to_atoms(i)
            writer.write(atoms = sys)

        return
    
    def save(self, fileName):
        np.save(fileName, self.gridStack)
        return

    def load(self, fileName):
        """Loads and overwrites the current grids as well as gridStack when given a valid filename (file extension not needed)"""
        # Load the gridStack from file
        self.gridStack = np.load(f"{fileName}.npy")

        # Then overwrite the current grids and simulation times with the loaded gridStack
        self.grid_S = self.gridStack[-1][0]
        self.grid_Mo = self.gridStack[-1][1]
        self.total_sim_time = self.gridStack[-1][2]

        # Then update all the other required stuff
        self.S_init = np.sum(self.grid_S)
        self.energy_cutoff_mean_S = self.get_energy_cutoff("S",0) # Calculate the energy cutoff for our atom type (S)
        self.a_S = self.a("S") # Calculate 'a'
        self.b_cutoff_mean_S = self.get_b_cutoff("S",0) # Calculate the b-cutoff for our atom type (S)
        self.rate_constant_S = self.get_rate_constant("S",0) # Calculate the rate constant for the system

        return

    def get_velocity_20C(self, atomSymb):
        """Given an atomic symbol, returns a velocity from a normal distribution"""
        # First get the mass (in kg)
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo

        # Define the width at 20C in m/s
        width = m.sqrt(124.1640143*self.k_B_si/m_n)

        # Using a normal distribution, pull a velocity
        velocity = np.random.normal(0,width)

        return velocity

    def get_transferred_energy(self, b, atomSymb, velocity):
        """Calculates and returns the transferred energy in eV, given the b-value in Å and velocity in m/s"""
        # Convert the b-value from Å to m, and electron energy from eV to J
        b *= 10**(-10)
        E_e = self.electronKin * 1.602176621*10**(-19)
        c = self.speed_of_light_si

        if atomSymb == "S":
            m_n = self.m_S
            Q = 16
        elif atomSymb == "Mo":
            m_n = self.m_Mo
            Q = 42

        # Use the b-value to calculate the reflection angle
        theta = 2*m.atan((self.coulomb_k_si*Q) / (self.m_e*b*self.electron_velocity**2))

        # Use the reflection angle as well as velocity to calculate the transferred energy
        nominator = 2*(E_e*(2*(c**2)*self.m_e + E_e) + m.sqrt(E_e*(2*(c**2)*self.m_e + E_e))*m_n*velocity*c)*(1 - m.cos(theta)) + (m_n*velocity*c)**2
        denominator = 2*m_n*self.speed_of_light_si**2
        E_T = nominator/denominator

        # Then convert to eV
        E_T *= 6.241509125*10**18

        ##### OLD #####
        ## First calculate the momentum
        #p_trans = (2*self.get_reduced_mass(atomSymb) * self.electron_velocity) / m.sqrt((b*10**(-10))**2 * self.a_S**2 + 1)
        ##### OLD #####

        ##### OLD #####
        ## Then use the momentum to calculate the energy
        #E_T = 6.241509125*10**18*p_trans**2/(2*m_n)
        ##### OLD #####

        return E_T

    def get_rate_constant(self, atomSymb):
        """Calculates and returns the rate constant of the whole system for ONE type of atom in 1/s"""
        # Determine how many atoms can get hit (In this case the amount of atoms in the bottom layer)
        numberS = np.sum(self.grid_S)

        # Determine the rate constant
        if atomSymb == "S":
            rate_constant = self.dose_rate * m.pi * self.b_cutoff_mean_S**2 * numberS
        else:
            rate_constant = self.dose_rate * m.pi * self.get_b_cutoff(atomSymb, 0)**2 * numberS # 1/s

        return rate_constant

    def update_rate_constant_S(self):
        """Updates the variable self.rate_constant_S"""
        self.rate_constant_S = self.get_rate_constant("S")

        return

    def get_relativistic_electron_velocity(self):
        """Calculates and returns the speed of the electrons, considering relativistic effects"""
        v_rela = self.speed_of_light_si*m.sqrt(1 - 1/(self.electronKin/(511024.6653) + 1)**2)
        return v_rela
    
    def get_electron_velocity(self):
        """Calculates and returns the classical electron velocity in m/s"""
        v_rest = m.sqrt(2*self.electronKin*1.602176621*10**(-19)/(self.m_e))

        return v_rest
    
    def get_relativistic_electron_mass(self):
        """Calculates the relativistic mass of the electron and returns it in kg"""
        v = self.get_relativistic_electron_velocity()
        m_r = (self.m_e) / m.sqrt(1 - (v/self.speed_of_light_si)**2)
        return m_r
    
    def get_reduced_mass(self, atomSymb):
        """Calculates and returns the reduced mass of the electron"""
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo
        
        reduced_mass = self.m_e * m_n / (self.m_e + m_n)
        return reduced_mass
    
    def get_displacement_cross_section(self, testMode = False):
        """Calculates and returns the current displacement cross-section for the bottom layer of the S-grid"""
        curS = np.sum(self.grid_S)
        scatSection = (self.S_init - curS) / (self.S_init * self.dose_rate * self.total_sim_time)
        
        # Convert from Å^2 to barn
        scatSection *= 10**8
        
        if testMode == True:
            return (self.S_init - curS)/self.S_init
        return scatSection

    def a(self, atomSymb):
        """Calculates the 'a' parameter, and returns it in 1/m"""
        # First get the variables in order
        v_0 = self.get_relativistic_electron_velocity()
        m_e = self.get_relativistic_electron_mass()
        #m_e = self.m_e

        if atomSymb == "S":
            m_n = self.m_S
            Q = 16
        elif atomSymb == "Mo":
            m_n = self.m_Mo
            Q = 42

        # Now calculate a
        a = (v_0**2 * m_e) / (self.coulomb_k_si * (1) * Q * (m_e/m_n + 1)**3)

        return a
    
    def get_p_cutoff(self, atomSymb):
        p = self.get_b_cutoff(atomSymb)**2
        return p

    def get_b_cutoff(self, atomSymb, velocity):
        """Calculates and returns the cutoff value for b in Å (angstrom)"""
        # Find the lowest TD value, as to find the b cutoff (as E_T ~ 1/b**2)
        if self.overwriteTD != False:
            TD_min = self.overwriteTD
        else:
            TD_min = self.TDlib["Td"].min()
        E_max = self.get_energy_cutoff(atomSymb, velocity)

        # Since we are limited by E_max, check whether this TD_Min is higher than E_Max
        if TD_min > E_max:
            warnings.warn(f"The electron energy is too low to damage the structure - choose a higher energy! Using 0.0005 Å as b_cutoff!")
            return 0.0005
        E = TD_min

        # Now get the mass of the atomic nucleus of the corresponding atom
        # The following should (for maximum compatibility) be some library but for now it's just some if-else statements
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo
        

        # Get the reduced mass of our electrons, as well as their velocity
        m_r = self.get_reduced_mass(atomSymb)
        m_rela = self.get_relativistic_electron_mass()
        v_0 = self.get_relativistic_electron_velocity()

        # Get 'a' (also known as kappa)
        a = self.a(atomSymb)

        # Calculate the cutoff value for b
        b_cutoff = (1/a) * m.sqrt(((2*m_rela*v_0)**2 / (E*1.602176621*10**(-19)*2*m_n)) - 1)
        #b_cutoff = (1/a) * m.sqrt(((2*m_r*v_0)**2 / (E*1.602176621*10**(-19)*2*m_n)) - 1)

        # Convert it to Å, increase it by 25%, and return it
        b_cutoff = 1.25 * b_cutoff * 10**(10)
    
        # Minimum limit
        if b_cutoff < 0.0005:
            b_cutoff = 0.0005

        return b_cutoff

    def get_energy_cutoff(self, atomSymb, velocity):
        """Calculates and returns the energy cutoff (maximum) in eV, given an atomic type and its velocity"""
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo

        # Convert eV to J
        E_e = self.electronKin * 1.602176621*10**(-19)

        # Calculate the electron rest energy
        E_0 = self.m_e * self.speed_of_light_si**2

        # Calculate the maximum energy transfer
        nominator = (2*m.sqrt(E_e*(E_e + 2*E_0)) + m_n*velocity*self.speed_of_light_si)**2
        denominator = 2*m_n*self.speed_of_light_si**2
        E_max = nominator/denominator

        #Then convert it to eV
        E_max *= 6.241509125*10**18

        return E_max

    def time_step(self):
        """Increases the total sim time as well as the current sim time using the total rate constant for the system"""
        # Calculate how long passed
        timePassed = -1 * m.log(uniform(0.0000000001, 1.0))/self.rate_constant_S # Change in seconds

        # Total time ran update
        self.total_sim_time += timePassed

        # This simulation cycle update
        self.current_sim_time += timePassed

        return

    def get_total_dose(self):
        self.total_dose = self.dose_rate * self.total_sim_time
        return self.total_dose

    def get_fingerPrint(self, layer, a1, a2):
        """Calculates and returns the fingerprint of the S atom on the given layer of the S grid, with the coordinates (layer, a1, a2)"""
        ### First check if the layer is valid, and if there is actually an atom at the given coordinates
        if layer not in [0, 1, 2]:
            raise Exception(f"The layer {layer} is not a valid layer. Choose either 0, 1 or 2.")
        if self.grid_S[layer][a1][a2] == False:
            return None

        ### Get the opposite layer
        if layer == 0:
            otherLayer = 2
        elif layer == 1:
            otherLayer = None
        elif layer == 2:
            otherLayer = 0

        ### Determine the amount of atoms in the same layer
        # First get which atoms should actually be looked into
        S_coords = []
        if a1 != self.squareSize - 1:
            S_coords.append((a1+1, a2))
            if a2 != 0:
                S_coords.append((a1+1, a2-1))
        if a2 != self.squareSize - 1:
            S_coords.append((a1, a2+1))
            if a1  != 0:
                S_coords.append((a1-1, a2+1))
        if a1 != 0:
            S_coords.append((a1-1, a2))
        if a2 != 0:
            S_coords.append((a1, a2-1))
        
        # Irregardless of which layer we are in, we need to count the neighbors as we would usually (yes even for middle-layer atoms)
        nS_NNs = 0
        for e in S_coords:
            nS_NNs += self.grid_S[layer][e[0],e[1]]

        # For layer 0 and 2 we need to look at the atom in the opposite layer, at the same coordinates
        if layer == 0 or layer == 2:
            nS_NNs += self.grid_S[otherLayer][a1][a2]


        ### Determine the amount of Mo neighbors, and their common neighbors with out S atom of interest
        # For the top, bottom as well as middle-layer
        # NOTE: INTERACTIONS WITH EDGE MIDDLE-LAYER ATOMS TURNED OUT TO BE RATHER SIMPLE (EXCEPT IN EDGE-CASES. FUCK EDGE CASES. THERE'S A 50% CHANCE OF THEM NOT COUNTING THE ATOMS IN THE SAME LAYER IN SOME STRUCTURES.)
        nMo = 0
        nS_list = []
        if a1 != 0 and a2 != self.squareSize-1 and self.grid_Mo[a1-1, a2]:
            nMo += self.grid_Mo[a1-1, a2]
            nS = 0
            if layer == 0 or layer == 2:
                nS += self.grid_S[otherLayer][a1,a2]
            nS += self.grid_S[layer][a1,a2+1]
            nS += self.grid_S[layer][a1-1,a2+1]
            nS_list.append(nS)

        if a1 != 0 and a2 != 0 and self.grid_Mo[a1-1, a2-1]:
            nMo += self.grid_Mo[a1-1, a2-1]
            nS = 0
            if layer == 0 or layer == 2:
                nS += self.grid_S[otherLayer][a1,a2]
            nS += self.grid_S[layer][a1,a2-1]
            nS += self.grid_S[layer][a1-1,a2]
            nS_list.append(nS)
            
        if a1 != self.squareSize-1 and a2 != 0 and self.grid_Mo[a1, a2-1]:
            nMo += self.grid_Mo[a1, a2-1]
            nS = 0
            if layer == 0 or layer == 2:
                nS += self.grid_S[otherLayer][a1,a2]
            nS += self.grid_S[layer][a1+1,a2]
            nS += self.grid_S[layer][a1+1,a2-1]
            nS_list.append(nS)

        ### Create the fingerprint
        nS_list.sort(reverse = True)
        while len(nS_list) < 4:
            nS_list.append(0)
        fingerPrint = [nMo] + [nS_NNs] + nS_list

        return tuple(fingerPrint)
    
    def get_TD(self, finger, situation, sample = True):
        """
        Using two indices, calculate the fingerprint of the atom and return its TD value, using the TD library. If no TD value exists for the given fingerprint, log this (add to the missingTDs list) and return False.
        """

        # Check for overwrite
        if self.overwriteTD != False:
            return self.overwriteTD
        
        # Define the TD library as local
        #TDlib = self.TDlib[self.TDlib["situation"] == situation]

        # Run through the pandas dataframe, and check if there are any corresponding value
        if len(self.TDlib[self.TDlib["situation"] == situation][self.TDlib[self.TDlib["situation"] == situation][self.fingerPrint] == finger]) == 0:
            # If there are none, add the fingerPrint to missingTDs (if it is not there already) and return None
            if (finger, situation) not in self.missingTDs:
                self.missingTDs.append((finger, situation))
            return None
        
        # If there are at least one corresponding TD value, either take the average of all the values and return the value, or sample one of the values and return it
        elif sample == True:
            return float(self.TDlib[self.TDlib["situation"] == situation][self.TDlib[self.TDlib["situation"] == situation][self.fingerPrint] == finger]["Td"].sample())
        else:
            return self.TDlib[self.TDlib["situation"] == situation][self.TDlib[self.fingerPrint] == finger].mean()["Td"]
    
    def get_missing_TDs(self):
        """
        Returns a list of all the fingerprints missing a TD value.
        """
        return self.missingTDs
    
    def clear_missing_TDs(self):
        """Clears the missing TD list."""
        self.missingTDs = []
        return

    def create_system(self, squareSize, structType):
        structType = structType.lower()

        if structType == "square 100":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)
            self.grid_S[1][:][:] = False # Set the middle-layer to be false
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False # Also set the sulfur at (0,0) to be false, as it is not actually there in a square structure

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

        elif structType == "square 50":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Set the middle-layer to be false, except for the edges
            self.grid_S[1][1:squareSize-1,1:squareSize-1] = False

            # Set the edges of the top and bottom layer to be false
            self.grid_S[0] = False
            self.grid_S[0][1:-1,1:-1] = True
            self.grid_S[-1] = False
            self.grid_S[-1][1:-1,1:-1] = True

            # Also set the sulfur at (0,0) to be false, as it is not actually there in a square structure
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)
        
        elif structType == "square mix-v1":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Set the middle-layer to be false, except for one edge
            self.grid_S[1][1:squareSize,0:squareSize] = False

            # Set the edges of the top and bottom layer to be false for one edge
            self.grid_S[0] = False
            self.grid_S[0][1:squareSize,0:squareSize] = True
            self.grid_S[-1] = False
            self.grid_S[-1][1:squareSize,0:squareSize] = True

            # Also set the sulfur at (0,0) to be false, as it is not actually there in a square structure
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

        elif structType == "square mix-v2":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Set the middle-layer to be false, except for two edges
            self.grid_S[1][1:squareSize,1:squareSize] = False

            # Set the edges of the top and bottom layer to be false for one edge
            self.grid_S[0] = False
            self.grid_S[0][1:squareSize,1:squareSize] = True
            self.grid_S[-1] = False
            self.grid_S[-1][1:squareSize,1:squareSize] = True

            # Also set the sulfur at (0,0) to be false, as it is not actually there in a square structure
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

        elif structType == "square mix-v3":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Set the middle-layer to be false, except for two edges
            self.grid_S[1][1:squareSize,0:squareSize-1] = False

            # Set the edges of the top and bottom layer to be false for one edge
            self.grid_S[0] = False
            self.grid_S[0][1:squareSize,0:squareSize-1] = True
            self.grid_S[-1] = False
            self.grid_S[-1][1:squareSize,0:squareSize-1] = True

            # Also set the sulfur at (0,0) to be false, as it is not actually there in a square structure
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

        elif structType == "square mix-v4":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Set the middle-layer to be false, except for two edges
            self.grid_S[1][1:squareSize,1:squareSize-1] = False

            # Set the edges of the top and bottom layer to be false for one edge
            self.grid_S[0] = False
            self.grid_S[0][1:squareSize,1:squareSize-1] = True
            self.grid_S[-1] = False
            self.grid_S[-1][1:squareSize,1:squareSize-1] = True

            # Also set the sulfur at (0,0) to be false, as it is not actually there in a square structure
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

        elif structType == "triangle 100 mo-edges":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Now put everything above a line to False
            for l in range(self.grid_S.shape[0]):
                for x in range(squareSize):
                    for y in range(squareSize):
                        if x + y >= squareSize + 1:
                            self.grid_S[l][x][y] = False

            # Then set the middle layer to false
            self.grid_S[1][:][:] = False

            # Also set the sulfur at (0,0) to be false, as it shouldn't be there
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

            # Then put everything above a line to False
            for x in range(squareSize - 1):
                for y in range(squareSize - 1):
                    if x + y >= squareSize - 1:
                        self.grid_Mo[x][y] = False
        
        elif structType == "triangle 50 mo-edges":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Now put everything above a line to False
            for l in range(self.grid_S.shape[0]):
                for x in range(squareSize):
                    for y in range(squareSize):
                        if x + y >= squareSize + 2:
                            self.grid_S[l][x][y] = False

            # Now set the middle layer to false
            self.grid_S[1][:][:] = False

            # Then set the outer part of the middle layer to True
            for x in range(squareSize):
                for y in range(squareSize):
                    if x + y == squareSize:
                        self.grid_S[1][x][y] = True

            self.grid_S[1][0,:] = True
            self.grid_S[1][:,0] = True

            # Now do the opposite for the upper and lower layer
            for l in range(self.grid_S.shape[0]):
                for x in range(squareSize):
                    for y in range(squareSize):
                        if x + y >= squareSize and (l == 0 or l == 2):
                            self.grid_S[l][x][y] = False

            for l in (0, 2):
                self.grid_S[l][0,:] = False
                self.grid_S[l][:,0] = False

            # Also set the sulfur at (0,0) to be false, as it shouldn't be there
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

            # Then put everything above a line to False
            for x in range(squareSize - 1):
                for y in range(squareSize - 1):
                    if x + y >= squareSize - 1:
                        self.grid_Mo[x][y] = False
            
        elif structType == "triangle 100 s-edges":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Now put everything above a line to False
            for l in range(self.grid_S.shape[0]):
                for x in range(squareSize):
                    for y in range(squareSize):
                        if x + y <= squareSize - 2:
                            self.grid_S[l][x][y] = False

            # Then set the middle layer to false
            self.grid_S[1][:][:] = False

            # Also set the sulfur at (0,0), (squareSize, 0) and (0, squareSize) to be false, as they shouldn't be there
            for i in range(self.grid_S.shape[0]):
                self.grid_S[i][0,0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

            # Then put everything above a line to False
            for x in range(squareSize-1):
                for y in range(squareSize-1):
                    if x + y <= squareSize - 3:
                        self.grid_Mo[x][y] = False
            
        elif structType == "triangle 50 s-edges":
            # Create S grid
            self.grid_S = np.ones((3, squareSize, squareSize), dtype = bool) # First comes layer, then x and then y. So: (l, x, y)

            # Now put everything above a line to False
            for l in range(self.grid_S.shape[0]):
                for x in range(squareSize):
                    for y in range(squareSize):
                        if x + y <= squareSize - 2:
                            self.grid_S[l][x][y] = False

            # Now set the middle layer to false
            self.grid_S[1][:][:] = False

            # Then set the outer part of the middle layer to True
            for x in range(squareSize):
                for y in range(squareSize):
                    if x + y == squareSize - 1:
                        self.grid_S[1][x][y] = True

            self.grid_S[1][-1,:] = True
            self.grid_S[1][:,-1] = True

            # Now do the opposite for the upper and lower layer
            for l in range(self.grid_S.shape[0]):
                for x in range(squareSize):
                    for y in range(squareSize):
                        if x + y == squareSize - 1 and (l == 0 or l == 2):
                            self.grid_S[l][x][y] = False

            for l in (0, 2):
                self.grid_S[l][-1,:] = False
                self.grid_S[l][:,-1] = False

            # Also set the sulfur at (0,0) to be false, as it shouldn't be there
            for l in range(self.grid_S.shape[0]):
                self.grid_S[l][0][0] = False

            # Create Mo grid
            self.grid_Mo = np.ones((squareSize-1, squareSize-1), dtype = bool) # This one only contains a single layer of Mo atoms, and as such it is simply (x, y)

            # Then put everything above a line to False
            for x in range(squareSize-1):
                for y in range(squareSize-1):
                    if x + y <= squareSize - 3:
                        self.grid_Mo[x][y] = False
            
        else:
            raise Exception("Structure type not recognized!")
        
        # Now create the grid to keep track of removed atoms
        self.grid_Removed = np.zeros((3, squareSize, squareSize), dtype = bool)

        return
