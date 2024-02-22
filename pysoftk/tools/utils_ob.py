from .utils_func import *
from openbabel import openbabel as ob
from openbabel import pybel as pb

import numpy as np
import os

def ff_ob_relaxation(mol, FF="MMFF94", iter_ff=100, ff_thr=1.0e-6):
    """Setting up an openbabel FF for 
       geometry optimization.
      
    Parameters
    ===========

    mol: OBabel.Mol
      An user-provided OpenBabel moelcule.

    FF: class.str
      Selected Force Field. Options are MMFF94, UFF, GAFF.

    iter: class.int
      Number of iterations to be used for the FF. 

    Returns
    ========

    mol_new: Mol.OBabel.Mol
      Open Babel Molecule with optmised Geometry.

    """

    ff = ob.OBForceField.FindForceField(str(FF))
    ff.Setup(mol.OBMol)
    ff.ConjugateGradients(int(iter_ff), float(ff_thr))

    ff.GetCoordinates(mol.OBMol)

    return mol

def rotor_opt(mol, FF="MMFF94", rot_steps=125):
    """Setting up an openbabel FF for rotor 
    optimization.
      
    Parameters
    ===========

    mol: OBabel.Mol
      An user-provided OpenBabel moelcule.

    FF: class.str
      Selected Force Field. Options are MMFF94, UFF, GAFF.

    rot_steps: class.int
      Number of iterations to be used for the rotational 
      optimization.

    Returns
    ========

    mol_new: Mol.OBabel.Mol
      Open Babel Molecule with optmised Geometry.

    """
    ff = ob.OBForceField.FindForceField(str(FF))
    ff.Setup(mol.OBMol)

    ff.FastRotorSearch()
    
    # N cycles, each with M forcefield ops
    ff.WeightedRotorSearch(int(rot_steps),int(np.ceil(rot_steps/2.0))) 
    ff.GetCoordinates(mol.OBMol)

    return mol 

def check_bond_order(file_name):
    """Function to check and correct the bonds
       of a given file. 

       Parameters
       ------------

       file_name: class.str
          An user-provided file containing a molecule
          in an accepted openbabel format.

       Returns
       ---------

        None:

         A new file containing checked bond and bond order
         molecular system.
    """
    base_name, file_extension = os.path.splitext(file_name)
    format_ext=file_extension.split(".")[-1]
    
    inputfile=pb.readfile(str(format_ext), str(file_name))
    mol=next(inputfile)

    mol.OBMol.DeleteHydrogens()
    mol.OBMol.ConnectTheDots() 
    mol.OBMol.PerceiveBondOrders()
    mol.OBMol.AddHydrogens()

    new_name="_".join([str(base_name),"new"])
    output=pb.Outputfile(str(format_ext),
                         ".".join([new_name,str(format_ext)]),
                         overwrite=True)
    output.write(mol)
    output.close()
