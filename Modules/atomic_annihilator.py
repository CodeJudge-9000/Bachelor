### INITIALIZE ###
from ase import Atoms
from ase.optimize import QuasiNewton
from gpaw import GPAW
from gpaw import Mixer, MixerSum, MixerDif
from gpaw.scf import Energy, Eigenstates, Density, Forces
import numpy as np

### CALCULATOR ###
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


def remove_atoms(system, atomIndex = [], atomRemoveIndex, relax = True, overwriteCalc = True):
    """
    Function to remove a list of atoms, or single atom. If relax = True as per default, the system will relax using its calculator after the atoms have been removed.
    
    Input:
     - system: Atoms system (ase)
     - atomIndex: one or multiple indices to update [int/list/np array]
     - atomRemoveIndex: one or multiple indices of atoms to remove [int/list/np array]
     - relax: whether to relax or not [bool] <OPTIONAL>
     - overwriteCalc: whether or not to overwrite the current calculator with the default one [bool] <OPTIONAL>
    
    Output:
     - Updated index/indices [int/np array]
    
    """
    # Update calculator if nescessary
    if overwriteCalc == True:
        system.calc = localcalc
    
    # Check input type
    if type(atomRemoveIndex) == list or type(atomRemoveIndex) == np.ndarray:
        listRemBool = True
        atomRemoveIndex = np.array(atomRemoveIndex)
    else:
        listRemBool = False
    
    if type(atomIndex) == list: #or type(atomRemoveIndex) == numpy.ndarray:
        listAtomBool = True
        atomIndex = np.array(atomIndex)
    else:
        listAtomBool = False
    
    # Never trust the user
    if listRemBool == True: # If a list of atoms need to be removed
        if listAtomBool == True:
            # If a list of indices need to be updated
            assert (len(np.intersect1d(atomIndex, atomRemoveIndex)) == 0), "Can't update atomic indices given, as they are part of the list to be removed!"
        else:
            # If a single index need to be updated
            assert (atomIndex not in atomRemoveIndex), "Can't update atomic index given, as it is part of the list to be removed!"
    
    else: # If a single atom need to be removed
        if listAtomBool == True:
            # If a list of indices need to be updated
            assert (atomRemoveIndex not in atomIndex), "Can't update atomic indices given, as the atom to be removed is part of the list to be updated!"
        else:
            # If a single index need to be updated
            assert (atomIndex != atomRemoveIndex), "Can't update atomic index given, as it is is to be removed!"
    
    if listRemBool == True:
        # Sort, then reverse the list
        atomRemoveIndex = np.sort(atomRemoveIndex)
        atomRemoveIndex = np.flip(atomRemoveIndex)
    
    # Remove the atoms to be removed
    if listRemBool == True:
        for i in atomRemoveIndex:
            del system[i]
    else:
        del system[atomRemoveIndex]
    
    # Update the atomic index/indices, and return it/them
    if listRemBool == True and listAtomBool == False:
        for i in atomRemoveIndex:
            if i < atomIndex:
                atomIndex -= 1
    elif listRemBool == False and listAtomBool == False:
        if atomRemoveIndex < atomIndex:
            atomIndex -= 1
    
    elif listRemBool == True and listAtomBool == True:
        for l in range(len(atomIndex)):
            for i in atomRemoveIndex:
                if i < atomIndex[l]:
                    atomIndex[l] -= 1
    elif listRemBool == False and listAtomBool == True:
        for l in range(len(atomIndex)):
            if atomRemoveIndex < atomIndex[l]:
                atomIndex[l] -= 1
    
    # Relax the system, assuming it already have a calculator
    if relax == True:
        opt = QuasiNewton(system)
        opt.run()
    
    return atomIndex

