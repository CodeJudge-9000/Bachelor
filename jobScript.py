### Initialize ###
import sys
import os
from ase.optimize import QuasiNewton
from ase import Atoms
from gpaw import GPAW
from gpaw.scf import Energy, Eigenstates, Density, Forces
from gpaw import Mixer, MixerSum, MixerDif
from ase import Atoms
from ase.io.trajectory import TrajectoryWriter
from ase.constraints import FixAtoms
from ase.io import read
from math import sqrt
import numpy as np
from ase.visualize import view
from ase.io import Trajectory
import matplotlib.pyplot as plt
sys.path.insert(0, '/zhome/06/4/136859/bachelor/mos/modules')
from bsubmission import *
from TD_Calculator import TD_relax_line
from atomic_annihilator import remove_atoms
from FingerPrints import SimpleFinger
import csv
from ase.parallel import paropen


################################### REPLACEABLES ###########################
sysString = SYSTEM_NAME
ArrayToRemove = THIS_ARRAY
extraString = ARRAY_NUMBERS
theIndex = THE_INDEX
theRadius = THE_RADIUS

############################################################################
# Reading system

sys = read(sysString)
copySys = sys.copy()

idNew = remove_atoms(system = copySys, atomIndex = theIndex, atomRemoveIndex = ArrayToRemove, relax = True)

fingerPrint = SimpleFinger(system = copySys, index = idNew)

Td = TD_relax_line(system = copySys, atomNo = idNew, atomNoToName = theIndex, stringExtra = extraString, subFolder = "defaultFolder", radius = theRadius, vaccuum = 7, endPoint = 0.0, stepsize = 0.25, smart = True, strNumbRand = False, delTraj = False)
    
# Changing to subfolder

os.chdir("./data")


# Writing a dataFile.
# Writing a dataFile.
d = pd.DataFrame({'finger': [fingerPrint], 'Td': [Td], 'original_id': [theIndex], \
     'startsystemname' : [sysString], 'removedList' : [ArrayToRemove]})

d.to_csv(f'id_{theIndex}_result_removal_{extraString}.csv')