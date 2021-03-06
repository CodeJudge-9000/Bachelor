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
from Nearest_Neighbors import *
import random
import time
from FingerPrints import *


def make_list_name(theList):
    theString = ""
    for i in theList:
        theString = theString + str(i)
    return theString


def make_data_point(theIndex, givensystemString,theRadius):
    """
    Inputs:
    theIndex -> index of the atom we want to move [int]
    givensystemString -> the string of the system eg. eq_struct_3x3.py [str]
    theRadius -> search radius for line relaxer [int]
    
    Output:
    Will make output from jobScript
    """
    
    sys = read(givensystemString)
    copySys = sys.copy()

    # Randomize
    #random.shuffle(allNNIndexes)
    # Nope not right now at least
    
    # Make list of lists
    RemovalLists = get_all_combs_tot(system = copySys, index = theIndex)

    ## 

    # Variables that we will change
    JobName = "placeHolder"
    coresAsked = 16

    # Unchanged variables
    scriptName = "cores.py"
    numHosts = 1
    email = "mathias_lendal@hotmail.dk"
    mem = "0.5GB"
    maxMem = "1.2GB"
    wallTime = "72:00"

    # In this section we just make a jobname, a scriptname, the number of hosts and so on for the file

    ## Script to copy name, is important, because we create a new file instead
    # so in case something goes wrong we don't ruin our original file
    scriptToCopyName = "jobScript.py"


    # Creating the object, that can do its wonders
    sub = bsubmissions(JobName, coresAsked, scriptName, numHosts, email, mem, maxMem, wallTime)

    # In this example i make the cores change in a for loop, where we can change something inside code, and submit, can be
    # done with basis as well.
    for removals in RemovalLists:
    # Here we make or open our dummy script
        repscript = open(scriptName, "w+")
        # Defining newVar

        # Changing jobname, not necessary but it is okay
        sub.jobName = f"Removals_{len(removals)}_Id{theIndex}"


        # From here everything is more automatic    
        #Making the sub file
        sub.make_bsub_file()

        # Opening original file and make the changeable variable
        fpy = open(scriptToCopyName, "r+")
        changeable = fpy.read()
        fpy.close()

        """ Replacing old variable hardcode with a new one, must change from what one wants to achive """
        changeable = changeable.replace(f'SYSTEM_NAME', f'\"{givensystemString}\"')
        changeable = changeable.replace(f'THIS_ARRAY', f'{removals}')
        listName = make_list_name(removals)
        changeable = changeable.replace(f'ARRAY_NUMBERS', f'\"{listName}\"')
        changeable = changeable.replace(f'THE_INDEX', f'{theIndex}')
        changeable = changeable.replace(f'THE_RADIUS', f'{theRadius}')

        # Deleting everything inside our dummy script
        repscript.truncate(0)
        repscript.close()

        # Open script, write our changed code and save it 
        repscript = open(scriptName, "w+")
        repscript.write(changeable)
        repscript.close()

        # Sending the job in
        sub.do_submission()

        # Delete subfile and wait such that we don't to things to quickly for system to understand
        sub.delete_bsub_file()
        time.sleep(5)
     
    """Special case, Here the 0 removal case"""
    repscript = open(scriptName, "w+")
    # Defining newVar

    # Changing jobname, not necessary but it is okay
    sub.jobName = f"Removals_Id{theIndex}"


    # From here everything is more automatic    
    #Making the sub file
    sub.make_bsub_file()

    # Opening original file and make the changeable variable
    fpy = open(scriptToCopyName, "r+")
    changeable = fpy.read()
    fpy.close()

    """ Replacing old variable hardcode with a new one, must change from what one wants to achive """
    changeable = changeable.replace(f'SYSTEM_NAME', f'\"{givensystemString}\"')
    changeable = changeable.replace(f'THIS_ARRAY', f'{[]}')
    listName = make_list_name(removals)
    changeable = changeable.replace(f'ARRAY_NUMBERS', f'\"none\"')
    changeable = changeable.replace(f'THE_INDEX', f'{theIndex}')
    changeable = changeable.replace(f'THE_RADIUS', f'{theRadius}')

    # Deleting everything inside our dummy script
    repscript.truncate(0)
    repscript.close()

    # Open script, write our changed code and save it 
    repscript = open(scriptName, "w+")
    repscript.write(changeable)
    repscript.close()

    # Sending the job in
    sub.do_submission()

    # Delete subfile and wait such that we don't to things to quickly for system to understand
    sub.delete_bsub_file()
    time.sleep(5)    

def make_uniquie_points(theIndex, givensystemString):
    theRadius = 5
    system6 = read(givensystemString)
    fingers, uniquerems = get_unqiue_removal_combi(system6,theIndex)
    noRem = [[]]
    rems = noRem + uniquerems

    """
    Inputs:
    theIndex -> index of the atom we want to move [int]
    givensystemString -> the string of the system eg. eq_struct_3x3.py [str]
    theRadius -> search radius for line relaxer [int]

    Output:
    Wil make output from jobScript
    """

    sys = read(givensystemString)
    copySys = sys.copy()


    ##

    # Variables that we will change
    JobName = "placeHolder"
    coresAsked = 16

    # Unchanged variables
    scriptName = "cores.py"
    numHosts = 1
    email = "mathias_lendal@hotmail.dk"
    mem = "0.5GB"
    maxMem = "1.2GB"
    wallTime = "72:00"

    # In this section we just make a jobname, a scriptname, the number of hosts and so on for the file

    ## Script to copy name, is important, because we create a new file instead
    # so in case something goes wrong we don't ruin our original file
    scriptToCopyName = "jobScript.py"


    # Creating the object, that can do its wonders
    sub = bsubmissions(JobName, coresAsked, scriptName, numHosts, email, mem, maxMem, wallTime)

    # In this example i make the cores change in a for loop, where we can change something inside code, and submit, can be
    # done with basis as well.
    for removals in rems:
    # Here we make or open our dummy script
        repscript = open(scriptName, "w+")

        # Changing jobname, not necessary but it is okay
        sub.jobName = f"Removals_{len(removals)}_Id{theIndex}"

        # From here everything is more automatic    
        #Making the sub file
        sub.make_bsub_file()

        # Opening original file and make the changeable variable
        fpy = open(scriptToCopyName, "r+")
        changeable = fpy.read()
        fpy.close()

        """ Replacing old variable hardcode with a new one, must change from what one wants to achive """
        changeable = changeable.replace(f'SYSTEM_NAME', f'\"{givensystemString}\"')
        changeable = changeable.replace(f'THIS_ARRAY', f'{removals}')
        listName = make_list_name(removals)
        changeable = changeable.replace(f'ARRAY_NUMBERS', f'\"{listName}\"')
        changeable = changeable.replace(f'THE_INDEX', f'{theIndex}')
        changeable = changeable.replace(f'THE_RADIUS', f'{theRadius}')

        # Deleting everything inside our dummy script
        repscript.truncate(0)
        time.sleep(1)
        repscript.close()
        time.sleep(1)
        # Open script, write our changed code and save it 
        repscript = open(scriptName, "w+")
        time.sleep(2)
        repscript.write(changeable)
        repscript.close()
        time.sleep(2)
        # Sending the job in
        sub.do_submission()
        time.sleep(2)
        # Delete subfile and wait such that we don't to things to quickly for system to understand
        sub.delete_bsub_file()
        time.sleep(2)

        