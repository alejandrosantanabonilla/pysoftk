from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom, rdMolTransforms
from rdkit.Chem import rdDistGeom as molDG

import numpy as np
from pysoftk.tools.utils_func import *
from pysoftk.tools.utils_rdkit import *

class Bd:
    """
    A class for creating branched polymers from RDKit molecules.

    **Examples**

    **Note**

    - Requires the RDKit package.
    - Functions `atom_neigh`, `count_plholder`, `swap_hyd`, and `no_swap` are imported from the `pysoftk.tools.utils_rdkit` module.

    """
    
    __slots__ = ['core', 'arm', 'atom']

    def __init__(self, core, arm, atom):
       """
       Initializes the Bd class.

       **Parameters**

       core: rdkit.Chem.rdchem.Mol
             The core molecule.
       arm: rdkit.Chem.rdchem.Mol
             The arm molecule.
       atom: str
             The placeholder atom for merging.

       """
       
       self.core = core
       self.arm = arm
       self.atom = atom

    def merge_arms(self, core, arm, atom):
       """
        Attaches user-defined molecules to a core molecule using a placeholder.

       Parameters
       ------------

        core: rdkit.Chem.rdchem.Mol
             The core molecule.

        arm: rdkit.Chem.rdchem.Mol
             The arm molecule.
        atom: str

             The placeholder atom.


        Returns
        --------

        mol: rdkit.Chem.rdchem.Mol
            The combined molecule.

       """
       outmol=Chem.CombineMols(core, arm)
       AllChem.EmbedMolecule(outmol, useRandomCoords=True)

       new_bond=atom_neigh(outmol, str(atom))
    
       br_1,c_1=tuple(new_bond[0])
       br_2,c_2=tuple(new_bond[-1])

       rwmol = Chem.RWMol(outmol)     
       rwmol.AddBond(c_1, c_2, Chem.BondType.SINGLE)

       rwmol.RemoveAtom(br_2)
       rwmol.RemoveAtom(br_1)

       mol=rwmol.GetMol()

       return mol

    
    def branched_polymer(self, force_field="MMFF94", iter_ff=100, swap_H=True):
       """
        Creates a branched polymer.

        Parameters
        ------------


        force_field: str, optional (default="MMFF94")
                    The force field to use (MMFF or UFF).
        
        iter_ff: int, optional (default=100)
                 The number of iterations for force field optimization.

        swap_H: bool, optional (default=True)
                Specifies whether to swap the placeholder atom with hydrogen.


        Returns
        --------

        newMol_H: rdkit.Chem.rdchem.Mol
                 The generated branched polymer.

        Raises
        -------

        ValueError: If an invalid force field is specified.

       """
       core=self.core
       arm=self.arm
       atom=self.atom
       
       res=self.merge_arms(core, arm, str(atom))
       AllChem.EmbedMolecule(res, useRandomCoords=True)
    
       nm_ph=count_plholder(res, str(atom))
    
       for _ in range(int(nm_ph)):
           res=self.merge_arms(res, arm, str(atom))

       valid_force_fields = ("MMFF", "UFF", "MMFF94")
       if force_field not in valid_force_fields:
            raise ValueError(f"Invalid force field: {force_field}. Valid options are: {valid_force_fields}")

       # Automatically change ff if necessary:
       if force_field == "MMFF94":
           force_field = "MMFF"  # Change to default MMFF94 for MMFF

       # Validate and convert iterations and steps to integers:
       try:
           iter_ff = int(iter_ff)
       except ValueError:
           raise ValueError("iter_ff must be an integer.")

       if swap_H:
         mol=swap_hyd(res, iter_ff, str(atom), force_field)

       if not swap_H:
         mol=no_swap(res, iter_ff, force_field)
            
       return mol


    

