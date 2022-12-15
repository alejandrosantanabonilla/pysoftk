from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom, rdMolTransforms
from rdkit.Chem import rdDistGeom as molDG

import numpy as np
from .utils import *
from pysoftk.linear_polymer.utils import *

class Bd:
    """
    A class for creating a branched polymer from given RDKit molecules.
     
    Examples
    ---------

    Note
    -----

    RDKit package must be installed.
    """
    
    __slots__ = ['core', 'arm', 'atom']

    def __init__(self, core, arm, atom):
       """
       Initialize this class.
          
       """
       
       self.core = core
       self.arm = arm
       self.atom = atom

    def merge_arms(self, core, arm, atom):
       """ Function to attach user defined molecules to a provided molecular core using a placeholder.

       Parameters
       ----------

       core : rdkit.Chem.rdchem.Mol
            RDKit Mol object

       arm : rdkit.Chem.rdchem.Mol
            RDKit Mol object
 
       atom : str
            The placeholder atom to combine the molecules and form a new monomer

       Return
       -------

       mol : rdkit.Chem.rdchem.Mol
             RDKit Mol object
       """
       
       outmol=Chem.CombineMols(core, arm)
       AllChem.EmbedMolecule(outmol)

       new_bond=seek_plhold(outmol, str(atom))
    
       # New way to flatten a list using sum function
       at=sum(new_bond, [])

       br_1,c_1=tuple(at[0])
       br_2,c_2=tuple(at[-1])

       rwmol = Chem.RWMol(outmol)     
       rwmol.AddBond(c_1, c_2, Chem.BondType.SINGLE)

       rwmol.RemoveAtom(br_2)
       rwmol.RemoveAtom(br_1)

       mol=rwmol.GetMol()

       return mol

    
    def branched_polymer(self, iter_ff=100, FF="MMFF", swap_H=True):
       """Function to create branched polymers

       Parameters
       -----------

       FF: str
            Selected FF between MMFF or UFF

       relax_iterations: int  
            Number of iterations used for relaxing a molecular object. 

        swap_H: bool
             Indicates if the user defined atomic place holder is changed to a Hydrogen atom or remain as the used species.   

       Return
       -------

       newMol_H : rdkit.Chem.rdchem.Mol
             RDKit Mol object
       """

       core=self.core
       arm=self.arm
       atom=self.atom
       
       res=self.merge_arms(core, arm, str(atom))
       AllChem.EmbedMolecule(res)
    
       nm_ph=count_plholder(res, str(atom))
    
       for _ in range(int(nm_ph)):
           res=self.merge_arms(res, arm, str(atom))

       if swap_H == "True":
         mol=swap_hyd(res, iter_ff, str(atom), FF)

       else:
         mol=no_swap(res, iter_ff, FF)
            
       return mol


    

