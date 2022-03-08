### INITIALIZE ###
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io.trajectory import TrajectoryWriter
from ase.optimize import QuasiNewton
from gpaw import GPAW
from gpaw import Mixer, MixerSum, MixerDif
from gpaw.scf import Energy, Eigenstates, Density, Forces
from math import sqrt
from ase.io import Trajectory
from pathlib import Path
import os
import numpy as np
import random



### DEFAULT CONVERGENCE CRITERIA ###
defaultConv = {'energy': 0.0005,  # eV / electron
     'density': 1.0e-4,  # electrons / electron
     'eigenstates': 4.0e-8,  # eV^2 / electron
     'bands': 'occupied'}

### CALCULATOR ###
#naming = "dist.txt"
localcalc = GPAW(mode='lcao',
            basis = 'szp(dzp)',
            xc='PBE',
            hund = False,
            #eigensolver = ETDM(searchdir_algo={'name': 'l-bfgs-p', 'memory': 10}),
            occupations={'name': 'fermi-dirac', 'width': 0.05},
            mixer=MixerDif(0.02, 5, 100),
            #nbands=16,
            symmetry='off' # We are moving stuff. Set this to 'off'
            #parallel={'domain': (5,2,2)}
            )

### Function to get distance between atoms ###
def atomDist(atom1, atom2):
    d1 = atom1.position[0] - atom2.position[0]
    d2 = atom1.position[1] - atom2.position[1]
    d3 = atom1.position[2] - atom2.position[2]
    dist = sqrt(d1**2 + d2**2 + d3**2)
    return dist

def distFromPoint(point, atom):
    d1 = point[0] - atom.position[0]
    d2 = point[1] - atom.position[1]
    d3 = point[2] - atom.position[2]
    dist = sqrt(d1**2 + d2**2 + d3**2)
    return dist

def TD_Line(system, ind, conv = defaultConv, step = 0.1, smart = True, UpDown = "up", endPoint = 0.0):
    """
    By pulling the atom out of its structure with a given step size, the TD is found as the energy required to remove the atom. The atom is pulled 10Å away from its structure on the third axis.
    
    By default, IF smart = False is set(!), the function will compute the TD by going in the positive direction, though if UpDown == "down" (default: "up") then it will compute by going in the negative direction. Do be aware that this will even flip the argument of endPoint to negative or positive accordingly, if given.
    
    By default the endPoint will be determined as 80% of the distance from the chosen atom to the wall of the cell in the z-direction. This can be overwritten by giving an argument different from endPoint = 0.0, but do be aware that the sign of the endpoint depends on whether smart is set to true or false, and then if smart == False it will depend on whether UpDown is "up" or "down".
    
    If smart == True (default), then UpDown will be ignored and the direction will be determined by whether the atom is above or below the plane for the center of mass. If above, go in positive z, if below go in negative z. If for some reason the atom is in the center of mass, it will default to positive.
    
    Input:
        - ase Atoms system [system]
        - int index for atom to calculate TD for [ind]
        - custom convergence criteria for LCAO [conv] <OPTIONAL>
        - stepsize in float Å [step] <OPTIONAL>
        - smart mode on/off [Bool] <OPTIONAL>
        - UpDown to determine positive or negative z-direction ("up" == positive, "down" == negative) [string] <OPTIONAL>
        - endPoint to determine how far the atom should be pulled, in Å [float] <OPTIONAL> (only enter a value if you're sure you know what you are doing)
    
    Output:
        - TD energy in eV
    
    """

    ### Initialize from 'smart' and 'UpDown' ###
    # Input sanitation
    step = abs(step)
    endPoint = abs(endPoint)
    UpDown = UpDown.lower()
    
    # Get center of mass
    com = system.get_center_of_mass()
    
    # Determine endpoint
    if endPoint != 0.0:
        endPoint = endPoint
    else:
        
        # distZ = (atom z - com z)
        # atom_cell_wall = cell_length/2 - distZ
        # endPoint = atom_cell_wall*0.8
        
        wallDist = system.get_cell()[-1][-1]/2 - (system[ind].position[2] - com[2])
        endPoint = wallDist*0.8
        
    
    if smart == True:
        step = abs(step)
        endPoint = abs(endPoint)
        
        # Determine if atom is above, below, or in the plane of the center of mass
        tempPos = system[ind].position
        if tempPos[2] > com[2] or tempPos[2] == com[2]:
            # If above or in the center of mass, leave the step and the endpoint as positive
            pass
        else:
            # Else (if it is then negative) go in the negative z-direction
            step *= -1
            endPoint *= -1
            
    elif UpDown == "down":
        # If the direction to calculate is "Down" in the z-direction, set the stepsize and endPoint to be negative
        step = abs(step)
        endPoint = abs(endPoint)
        step *= -1
        endPoint *= -1
            
    elif UpDown == "up":
        # If the direction to calculate is "up" in the z-direction, set the stepsize and endPoint to be positive (even though the stepsize is already positive, this is to avoid bugs by future modification to the code)
            step = abs(step)
            endPoint = abs(endPoint)
    
    
    ### Initialize systems & get reference energy ###
    # Initialize system to be 'damaged'
    sysCopy = system.copy()
    sysCopy.calc = localcalc  # Vi skal passe på her, fordi der er valgt en localcalc vi ikke må vælge om, uden at modifie modulet.
    
    # Reference energy
    system.calc = localcalc
    eRef = system.get_total_energy()



    ### Calculate and write relative energy for different distances ###

    # Arrays to store data in
    eDifs = np.array([])
    dists = np.array([])

    # Remove/Damage the system material by pulling away the atom by the displacement until reaching 10Å
    for disp in np.arange(0,endPoint+step,step):
        sysCopy[ind].position[2] += step # The third coordinate is the z-coordinate, and {step}Å is added for each iteration


        # Now calculate the following energy difference: Damaged System - Undamaged(reference) system
        eCopy = sysCopy.get_total_energy()
        eDif = eCopy - eRef

        # Update arrays
        eDifs = np.append(eDifs, eDif)
        dists = np.append(dists, disp)
        
    #return the largest energy difference from zero (eRef) as TD
    return np.max(eDifs)




def TD_Quick(system, ind, conv = defaultConv, radius = 5, smart = True, UpDown = "up", endPoint = 0.0):
    """
    Alternative version of TD_line, where only the first and last point at 10Å is calculated. As with TD_Line, there is the UpDown and smart argument to determine the direction to calculate TD.
    
    By default, IF smart = False is set(!), the function will compute the TD by going in the positive direction, though if UpDown == "down" (default: "up") then it will compute by going in the negative direction. Do be aware that this will even flip the argument of endPoint to negative or positive accordingly, if given.
    
    By default the endPoint will be determined as 80% of the distance from the chosen atom to the wall of the cell in the z-direction. This can be overwritten by giving an argument different from endPoint = 0.0, but do be aware that the sign of the endpoint depends on whether smart is set to true or false, and then if smart == False it will depend on whether UpDown is "up" or "down".
    
    If smart == True (default), then UpDown will be ignored and the direction will be determined by whether the atom is above or below the plane for the center of mass. If above, go in positive z, if below go in negative z. If for some reason the atom is in the center of mass, it will default to positive.
    
    
    Input:
    - ase Atoms system [system]
    - int index for atom to calculate TD for [ind]
    - custom convergence criteria for LCAO [conv] <OPTIONAL>
    - radius for sphere within NOT to clamp atoms [number] <OPTIONAL>
    - smart mode on/off [Bool] <OPTIONAL>
    - UpDown to determine positive or negative z-direction ("up" == positive, "down" == negative) [string] <OPTIONAL>
    - endPoint to determine how far the atom should be pulled, in Å [float] <OPTIONAL> (only enter a value if you're sure you know what you are doing)
    
    Output:
        - TD energy in eV
    """
    ### Initialize from 'smart' and 'UpDown' ###
    # Input sanitation
    endPoint = abs(endPoint)
    UpDown = UpDown.lower()
    
    # Get center of mass
    com = system.get_center_of_mass()
    
    # Determine endpoint
    if endPoint != 0.0:
        endPoint = endPoint
    else:
        
        # distZ = (atom z - com z)
        # atom_cell_wall = cell_length/2 - distZ
        # endPoint = atom_cell_wall*0.8
        
        wallDist = system.get_cell()[-1][-1]/2 - (system[ind].position[2] - com[2])
        endPoint = wallDist*0.8
        
    
    if smart == True:
        endPoint = abs(endPoint)
        
        # Determine if atom is above, below, or in the plane of the center of mass
        tempPos = system[ind].position
        if tempPos[2] > com[2] or tempPos[2] == com[2]:
            # If above or in the center of mass, leave the endpoint as positive
            pass
        else:
            # Else (if it is then negative) go in the negative z-direction
            endPoint *= -1
            
    elif UpDown == "down":
        # If the direction to calculate is "Down" in the z-direction, set the endPoint to be negative
        endPoint = abs(endPoint)
        endPoint *= -1
            
    elif UpDown == "up":
        # If the direction to calculate is "up" in the z-direction, set the endPoint to be positive
            endPoint = abs(endPoint)
    
    if smart == True:
        # Get center of mass
        com = system.get_center_of_mass()
        
        # Determine if atom is above, below, or in the plane of the center of mass
        tempPos = system[ind].position
        if tempPos[2] > com[2] or tempPos[2] == com[2]:
            # If above or in, do nothing
            pass
        else:
            # Else go in the negative z-direction
            endPoint *= -1
            
    elif UpDown == "down":
        # If the direction to calculate is "Down" in the z-direction, set the endPoint to be negative
            endPoint *= -1
    
    
    ### Initialize systems & get reference energy ###
    # Initialize system to be 'damaged'
    sysCopy = system.copy()
    sysCopy.calc = localcalc
    opt = QuasiNewton(sysCopy)
    
    #temp
    tw = TrajectoryWriter(f"TDtest_{sysCopy[ind].symbol}.traj", mode='w', atoms=sysCopy)
    tw.write()
    
    # Reference energy
    system.calc = localcalc
    eRef = system.get_total_energy()

    
    ### Move atom to endPoint, and relax system with constriction around initial point ###
    constraints = FixAtoms(indices=[atom.index for atom in sysCopy \
                                            if (atomDist(sysCopy[ind],atom) >= radius or
                                                atomDist(sysCopy[ind],atom) <= 0.001)])
    sysCopy.set_constraint(constraints)
    sysCopy[ind].position[2] += endPoint
    opt.run()
    sysCopy.set_constraint()
    tw.write()
    
    ### Calculate and write relative energy for endpoint, with  ###
    eCopy = sysCopy.get_total_energy()
    eDif = eCopy - eRef
    
    return eDif




def TD_relax_line(system, atomNo, atomNoToName = None, stringExtra = "", subFolder = "defaultFolder", radius = 5, vaccuum = 7, endPoint = 0.0, stepsize = 0.25, smart = True, strNumbRand = True, delTraj = False):
    """
    Input:
    - system: atoms object
    - atomNo: what atom number should be pulled
    - atomNoToName: number of atom to add to the name of the .traj file
    - stringExtra: String which is added to the name of the trajectory file
    - radius: The radius from the atom which should not be clamped
    - vaccuum: The amount of vacuum we want to put
    - endPoint: Distance from
    - stepsize: Stepsize for the line relaxation
    - smart: bool indicating whether the direction of pulling should be determined automatically (affects stepsize and endPoint)
    - strNumbRand: bool indicating whether a random number should be added to the trajectory file name (important in case multiple runs are done on the same atom)
    - keepTraj: bool indicating whether or not the trajectory file should be deleted after use
    
    Output:
    - Trajectory file in folder
    - TD value

    """
    
    ### Initialize from 'smart' & vaccuum ###
    # Centering system
    vaccuum = abs(vaccuum)
    system.center(vacuum = vaccuum)
    
    # Get center of mass
    com = system.get_center_of_mass()

    # Determine endpoint
    if endPoint != 0.0:
        endPoint = endPoint
    else:
        # distZ = (atom z - com z)
        # atom_cell_wall = cell_length/2 - distZ
        # endPoint = atom_cell_wall*1
        
        wallDist = system.get_cell()[-1][-1]/2 - (system[atomNo].position[2] - com[2])
        endPoint = wallDist*1
    
    if smart == True:
        stepsize = abs(stepsize)
        endPoint = abs(endPoint)
        
        # Determine if atom is above, below, or in the plane of the center of mass
        tempPos = system[atomNo].position
        if tempPos[2] >= com[2]:
            # If above or in the center of mass, leave the stepsize and the endpoint as positive
            pass
        else:
            # Else (if it is then negative) go in the negative z-direction
            stepsize *= -1
            endPoint *= -1
    
    ############################# OLD FOR STEPSIZE AS LIST #########################
    #if smart == True:
    #    stepsizes = [abs(e) for e in stepsizes]
    #    endPoint = abs(endPoint)
    #    
    #    # Determine if atom is above, below, or in the plane of the center of mass
    #    tempPos = system[atomNo].position
    #    if tempPos[2] > com[2] or tempPos[2] == com[2]:
    #        # If above or in the center of mass, leave the stepsizes as positive
    #        pass
    #    else:
    #        # Else (if it is then negative) go in the negative z-direction
    #        stepsizes = [-1*e for e in stepsizes]
    #        endPoint *= -1
    ################################################################################

    ### Pulling and Relaxation ###
    id = atomNo
    origPos = system[id].position[2]
    
    # This is a secret list we'll use for later ;)
    origIndices=[atom.index for atom in system \
                                    if (atomDist(system[id],atom) <= radius and
                                        atomDist(system[id],atom) >= 0.001)]

    # Make copy of the system #
    systemCopy = system.copy()
    systemCopy.calc = localcalc # Setting system to local calculator
    opt = QuasiNewton(systemCopy)
    
    # If we need to keep the .traj files, create a subfolder if it does not exist
    if delTraj == False:
        path = os.getcwd()
        Path(f"{path}/{subFolder}").mkdir(parents=True, exist_ok=True)
        #os.system(f"mkdir {subFolder}")
        
        # Also change the working directory to the folder
        os.chdir(f"{path}/{subFolder}")

    # Look at atomNoToName and update it if nescessary
    if atomNoToName == None:
        atomNoToName = atomNo
    
    # Check if a random number should be added
    if strNumbRand == True:
        randString = f"_{random.randrange(0, 100000)}"
    else:
        randString = ""
        
    # Check if an extra string should be added, and modify it if true
    if stringExtra != "":
        stringExtra = f"_{stringExtra}"
    
    # Open file
    string = f"pull_relax_r{radius}_s{abs(stepsize)}_atomNo{atomNoToName}_{system[atomNo].symbol}{stringExtra}{randString}.traj"
    tw = TrajectoryWriter(string, mode='w', atoms=systemCopy)


    # Get the initial system energy/reference energy, and create a list to add to
    energies = []
    E = systemCopy.get_total_energy()
    energies.append(E)
    
    for c in np.arange(0, endPoint+stepsize, stepsize):
        # Pull
        systemCopy[id].position[2] = origPos + c

        # Save systemCopy
        #tw.write()  Removed for now to get more easy relax energies

        # Restrict atoms & use the secret list
        indices=[atom.index for atom in systemCopy \
                                        if (atomDist(systemCopy[id],atom) >= radius or
                                            atomDist(systemCopy[id],atom) <= 0.001)]

        for e in origIndices:
            if e in indices:
                indices.remove(e)

        constraints = FixAtoms(indices)
        systemCopy.set_constraint(constraints)

        # Relax
        opt = QuasiNewton(systemCopy)
        opt.run(fmax = 0.05)

        # Unrestrict systemCopy
        systemCopy.set_constraint()
        
        # Get the system energy, and append it
        E = systemCopy.get_total_energy()
        energies.append(E)

        # Save systemCopy
        tw.write()
    
    # Find the TD value, and return it
    energies = np.array(energies)
    TD = np.max(energies) - energies[0]
    
    # Delete trajectory file, if needed
    if delTraj == True and os.path.exists(string):
        os.remove(string)
    # Otherwise, return to the original working directorys
    elif delTraj == False:
        os.chdir(path)
    
    return TD



def TD_From_Traj(trajString):
    # Load system from .traj file
    readSys = Trajectory(trajString)
    
    # Create array of energies
    arr = np.array([])
    for i in range(len(readSys)):
        arr = np.append(arr, readSys[i].get_total_energy())

    # Calculate TD
    TD = np.max(arr) - arr[0] # TD as the highest energy minus the initial (reference) energy
    
    return TD




