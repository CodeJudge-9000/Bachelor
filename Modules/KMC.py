
# Initialization
from ase import Atoms
from ase.io.trajectory import TrajectoryWriter
import warnings
import numpy as np
import math as m
from random import uniform
from random import randint


class KMC:
    def __init__(self, gridSize, structureType, TDlib, fingerPrint, kinetic_E, electron_dose):
        # First create the system
        self.create_system(gridSize, structureType)

        # Define a bunch of internal parameters
        self.squareSize = gridSize
        self.structType = structureType
        self.dose = electron_dose # number of electrons/(Å**2 * s)
        self.total_sim_time = 0 # Some time unit. Figure this one out later
        self.fingerPrint = fingerPrint
        self.TDlib = TDlib # TD energies in eV
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV) <- Units are *important*
        self.gridStack = np.array([np.array([self.grid_S, self.grid_Mo, self.total_sim_time], dtype=object)]) # np array used to store the two grids, as well as the 
        self.S_init = np.sum(self.grid_S[-1])

        # Define some constants
        self.m_e = 9.1093837015*10**(-31) # Electron mass in kg
        self.m_S = 32.06 * 1.660540200*10**(-27) # Sulfur mass in kg
        self.m_Mo = 95.95 * 1.660540200*10**(-27) # Molybdenum mass in kg
        self.speed_of_light_si = 2.998*10**8 # 'c' in m/s
        self.coulomb_k_si = 2.307077515*10**(-28) # The coulomb constant times e^2, given in SI units (kg*m^3)/(s^2)

        # Calculate some of the constant values of the system
        self.rate_constant_S = self.get_rate_constant("S") # Calculate the rate constant for the system
        self.relativistic_electron_mass = self.get_relativistic_electron_mass() # The relativistic electron mass in kg
        self.electron_velocity = self.get_electron_velocity() # Electron velocity in m/s
        self.energy_cutoff_S = self.get_energy_cutoff("S") # Calculate the energy cutoff for our atom type (S)
        self.a_S = self.a("S") # Calculate 'a'
        self.b_cutoff_S = self.get_b_cutoff("S") # Calculate the b-cutoff for our atom type (S)

        # Create the initial (and empty) missing fingerprints list
        self.missingTDs = [] # Add the missing fingerprints to this list
    


    def set_electron_energy(self, kinetic_E):
        """Updates the electron energy"""
        self.electronKin = kinetic_E # Kinetic energy of electrons. Same unit as TD values (eV)

        # Remember to update the dependent functions
        self.electron_velocity = self.get_electron_velocity()
        self.energy_cutoff_S = self.get_energy_cutoff("S")
        self.relativistic_electron_mass = self.get_relativistic_electron_mass()
        self.a_S = self.a("S")
        self.b_cutoff_S = self.get_b_cutoff("S")
        return
    
    def set_electron_dose(self, electron_dose):
        """Updates the electron dose"""
        self.dose = electron_dose

        # Remember to update the dependent functions
        self.rate_constant_S = self.get_rate_constant("S")
        return

    def run(self, runTime):
        self.current_sim_time = 0
        
        while self.current_sim_time < runTime:
            self.simulate_electron()
            self.time_step()
            self.gridStack = np.concatenate((self.gridStack, np.array([np.array([self.grid_S, self.grid_Mo, 1],dtype=object)])))
        
        return
    
    def simulate_electron(self): # Fix this thing
        # Choose which electron to interact with
        sideLen = len(self.grid_S[-1])
        a1, a2 = randint(0, sideLen-1), randint(0, sideLen-1)

        # Get the fingerprint for the atom
        # First in the bottom layer
        fingerPrint = self.get_fingerPrint(2, a1, a2)
        if fingerPrint == None:

            # Then if there is no atom in the bottom layer, check if there is one in the middle layer
            fingerPrint = self.get_fingerPrint(1, a1, a2)
            if fingerPrint == None:
                return 1

        # Figure out the interaction distance (b-value) of the electron and atom
        b = m.sqrt(uniform(0, self.b_cutoff_S**2))
        
        # Figure out how much energy is transferred
        E_T = self.get_transferred_energy(b, "S")

        #print(f"E_T: {E_T}")
        #print(f"self.energy_cutoff_S: {self.energy_cutoff_S}")
        if E_T > self.energy_cutoff_S:
            E_T = self.energy_cutoff_S
        
        # Now check whether the transferred energy is higher than the TD value for this atom
        TD = self.get_TD(fingerPrint)

        #print(f"TD: {TD}")
        if TD == None:
            # Return 1 if there is no corresponding TD value
            return 1
        
        if E_T >= TD:
            # remove atom if the transferred energy exceeds the TD value
            self.grid_S[-1][a1][a2] = False

            # Then update the rate constant for the system
            self.rate_constant_S = self.get_rate_constant()

            # Then return 0
            return 0
        else:
            # Otherwise, return 1
            return 1

        return 1
    
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

        return
    
    def reset(self):
        """Creates and updates the current grids with undamaged ones using the size of the current grids. Also resets the total simulation time."""
        self.create_system(self.squareSize, self.structType)
        self.total_sim_time = 0

        return

    def current_grid_to_atoms(self):
        """Converts the current grids into an ase atoms object, and returns it"""
        return self.grid_to_atoms(-1)

    def grid_to_atoms(self, stackLayer):
        """Converts the given gridStack grids into an ase atoms object, and returns it"""
        # Get the S and Mo grid
        grid_S = self.gridStack[stackLayer][0]
        grid_Mo = self.gridStack[stackLayer][1]

        atomS_list = []
        atomMo_list = []
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
                    if grid_S[L][a1][a2] == True:
                        atomS_list.append([a1*intera1 + a2*intera2[0],a2*intera2[1],(L-1)*interL])

        # Now construct the list of Mo atoms, using much the same method
        # Go over each first coordinate (aka row 'r' or a1)
        for a1 in range(grid_Mo.shape[0]):
            # Go over each second coordinate (aka column 'c' or a2)
            for a2 in range(grid_Mo.shape[1]):
                # Begin at the first column in the first layer
                if grid_Mo[a1][a2] == True:
                    atomMo_list.append([a1*intera1 + a2*intera2[0] + 3.18,a2*intera2[1] + 1.836, 0])

        # Now that we've constructed the coordinate lists, we construct the name of the system
        sysString = f"S{len(atomS_list)}Mo{len(atomMo_list)}"

        # And now we finally construct the Atoms object
        system = Atoms(sysString,
                    positions=atomS_list + atomMo_list,
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
        self.S_init = np.sum(self.grid_S[0])
        self.rate_constant_S = self.get_rate_constant("S") # Calculate the rate constant for the system
        self.energy_cutoff_S = self.get_energy_cutoff("S") # Calculate the energy cutoff for our atom type (S)
        self.a_S = self.a("S") # Calculate 'a'
        self.b_cutoff_S = self.get_b_cutoff("S") # Calculate the b-cutoff for our atom type (S)

        return

    def get_transferred_energy(self, b, atomSymb):
        """Calculates and returns the transferred energy in eV, given the b-value in Å"""
        # First calculate the momentum
        p_trans = (2*self.get_reduced_mass(atomSymb) * self.electron_velocity) / m.sqrt((b*10**(-10))**2 * self.a_S**2 + 1)

        # Then find other required parameters
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo

        # Then use the momentum to calculate the energy
        E_T = 6.241509125*10**18*p_trans**2/(2*m_n)

        return E_T

    def get_rate_constant(self, atomSymb):
        """Calculates and returns the rate constant of the whole system for ONE type of atom in 1/s"""
        # Determine how many atoms can get hit (In this case the amount of atoms in the bottom layer)
        numberS = np.sum(self.grid_S[-1])

        # Determine the rate constant
        rate_constant = self.dose * m.pi * self.get_b_cutoff(atomSymb)**2 * numberS # 1/s
        #print("b_cutoff: ",self.get_b_cutoff(atomSymb))
        #print("rate_constant: ", rate_constant)

        return rate_constant

    def get_relativistic_electron_velocity(self):
        """Calculates and returns the speed of the electrons, considering relativistic effects"""
        v_rela = 2.998*10**8*m.sqrt(1 - 1/(self.electronKin/(5.1098895*10**5) + 1)**2)
        return v_rela
    
    def get_electron_velocity(self):
        """Calculates and returns the classical electron velocity in m/s"""
        v_rest = m.sqrt(2*self.electronKin*1.602176621*10**(-19)/(self.m_e))

        return v_rest
    
    def get_relativistic_electron_mass(self):
        """Calculates the relativistic mass of the electron and returns it in kg"""
        v = self.get_electron_velocity()
        m_r = (self.m_e) / (1 - (v/self.speed_of_light_si)**2)
        return m_r
    
    def get_reduced_mass(self, atomSymb):
        """Calculates and returns the reduced mass of the system"""
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo
        
        reduced_mass = self.m_e * m_n / (self.m_e + m_n)
        return reduced_mass
    
    def get_scattering_cross_section(self):
        """Calculates and returns the current scattering cross-section for the bottom layer of the S-grid"""
        curS = np.sum(self.grid_S[-1])
        scatSection = (self.S_init - curS) / (self.S_init * self.dose * self.total_sim_time)

        return scatSection

    def a(self, atomSymb):
        """Calculates the 'a' parameter, and returns it in 1/m"""
        # First get the variables in order
        v_0 = self.get_electron_velocity()

        if atomSymb == "S":
            m_n = self.m_S
            Q = 8
        elif atomSymb == "Mo":
            m_n = self.m_Mo
            Q = 21

        # Now calculate a
        a = (v_0**2 * self.m_e) / (self.coulomb_k_si * (1) * Q * (self.m_e/m_n + 1)**3)

        return a
    
    def get_p_cutoff(self, atomSymb):
        p = self.get_b_cutoff(atomSymb)**2
        return p

    def get_b_cutoff(self, atomSymb):
        """Calculates and returns the cutoff value for b in Å (angstrom)"""
        # Find the lowest TD value, as to find the b cutoff (as E_T ~ 1/b**2)
        TD_min = self.TDlib["Td"].min() * 1.05
        E_max = self.get_energy_cutoff(atomSymb)

        # Since we are limited by E_max, check whether this TD_Min is higher than E_Max
        if TD_min > E_max:
            E = E_max
            warnings.warn(f"The electron energy is too low to damage the structure - choose a higher energy!")
        else:
            E = TD_min

        # Now get the mass of the atomic nucleus of the corresponding atom
        # The following should (for maximum compatibility) by some library but for now it's just some if-else statements
        if atomSymb == "S":
            m_n = self.m_S
        elif atomSymb == "Mo":
            m_n = self.m_Mo
        

        # Get the relativistic mass of our electrons, as well as their velocity
        m_r = self.get_reduced_mass(atomSymb)
        v_0 = self.get_relativistic_electron_velocity()

        # Get 'a'
        a = self.a(atomSymb)

        # Calculate the cutoff value for b
        b_cutoff = (1/a) * m.sqrt((2*m_r*v_0)**2 / (E*1.602176621*10**(-19)*2*m_n*1.660540200*10**(-27)) - 1)

        # Convert it to Å and return it
        b_cutoff = b_cutoff * 10**(10)
        #print(b_cutoff)

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
        timePassed = -1 * m.log(uniform(0.0000000001, 1.0))/self.rate_constant_S # Change in seconds

        # Total time ran update
        self.total_sim_time += timePassed

        # This simulation cycle update
        self.current_sim_time += timePassed

        return

    def get_fingerPrint(self, layer, a1, a2):
        """Calculates and returns the fingerprint of the S atom on the given layer of the S grid, with the coordinates (layer, a1, a2)"""
        # First check if the layer is valid, and if there is actually an atom at the given coordinates
        if layer not in [0, 1, 2]:
            raise Exception(f"The layer {layer} is not a valid layer. Choose either 0, 1 or 2.")
        if self.grid_S[layer][a1][a2] == False:
            return None

        # Get the opposite layer
        if layer == 0:
            otherLayer = 2
        elif layer == 1:
            otherLayer = None
        elif layer == 2:
            otherLayer = 0

        # Get the amount of NN S atoms there are in the same layer
        nS_NNs = 0
        for e in [(a1+1, a2), (a1+1, a2-1), (a1, a2-1), (a1, a2+1), (a1-1, a2+1), (a1-1, a2)]:
            nS_NNs += self.grid_S[layer][e[0],e[1]]

        # If we are in the middle layer (1) then we will need to count in the two other layers (using the same method as above)
        if layer == 1:
            for e in [(a1+1, a2), (a1+1, a2-1), (a1, a2-1), (a1, a2+1), (a1-1, a2+1), (a1-1, a2)]:
                nS_NNs += self.grid_S[0][e[0],e[1]]
                nS_NNs += self.grid_S[2][e[0],e[1]]

        # Otherwise, we need to look at the atom in the opposite layer, at the same coordinates
        if layer == 0 or layer == 2:
            nS_NNs += self.grid_S[otherLayer][a1][a2]

        # Now determine the amount of Mo neighbors (and save the coordinates in a list)
        # Given the coordinates of a S atom, the Mo neighbor would be at:
        #   (a1 - 1, a2 - 1), (a1 - 1, a2), (a1, a2 - 1)
        toReview = [(a1 - 1, a2 - 1, 0)] # a1 and a2 will never be 0 at the same time for our square molecule, and the corner at max coords is no issue
        if a1 != 0 and a2 != len(self.grid_Mo):
            toReview.append((a1 - 1, a2, 1)) # The second would be to the left of the first one
        if a2 != 0 and a1 != len(self.grid_Mo):
            toReview.append((a1, a2 - 1, 2)) # The third would be downwards from the first one
        
        # Now check how many Mo atoms are at the positions, and how many common neighbors there are
        nMo = 0
        nS_list = []
        #print("toReview: ",toReview)
        for e in toReview:
            nS = 0
            if self.grid_Mo[e[0],e[1]]:
                nMo += 1
                
                # For the middle layer
                if layer == 1:
                    raise NotImplementedError("Middle-layer atoms have yet to be implemented.")
                
                # For the top and bottom layer
                else:
                    nS += self.grid_S[layer][e[0]+1][e[1]]
                    nS += self.grid_S[layer][e[0]][e[1]+1]
                    nS += self.grid_S[otherLayer][e[0]+1][e[1]+1] # This is the same position as our atom, but in the opposite layer
            
            nS_list.append(nS)
        nS_list.sort(reverse = True)
        
        while len(nS_list) < 4:
            nS_list.append(0) # The fingerprint allows for up to 4 NN Mo atoms, but in our case there will at most be 3, so this is a quick fix

        # Now create the actual fingerprint
        fingerPrint = [nMo] + [nS_NNs] + nS_list

        return fingerPrint
    
    def get_TD(self, finger):
        """
        Using two indices, calculate the fingerprint of the atom and return its TD value, using the TD library. If no TD value exists for the given fingerprint, log this (add to the missingTDs list) and return False.
        """
        # First get the fingerprint for the corresponding indices
        #finger = self.get_fingerPrint(a1, a2)

        # Now run through the pandas dataframe, and check if there are any corresponding value
        if len(self.TDlib[self.TDlib[self.fingerPrint] == str(finger)]) == 0:
            # If there are none, add the fingerPrint to missingTDs (if it is not there already) and return True
            if finger not in self.missingTDs:
                self.missingTDs.append(finger)
            return None
        
        # If there are at least one corresponding TD value, take the average of all the values and return the value
        else:
            return self.TDlib[self.TDlib[self.fingerPrint] == str(finger)].mean()["Td"]
    
    def get_missing_TDs(self):
        """
        Returns a list of all the fingerprints missing a TD value.
        """
        return self.missingTDs
    
    def clear_missing_TDs(self):
        """Clears the missing TD list."""
        self.missingTDs = []
        return

