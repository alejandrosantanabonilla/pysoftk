import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import copy

from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer import super_monomer  as sm
from pysoftk.linear_polymer import mol_conformer  as mc
from pysoftk.linear_polymer.utils import *
from .utils import * 

class Rn:
    """A class for creating a circular (ring-shaped) polymer from given RDKit molecules.

    Examples
    ---------

    Note
    -----

    RDKit package must be installed.
    """
    
    __slots__ = ['mol', 'atom']

    def __init__(self, mol, atom):
       """Initialize this class.
       """
       
       self.mol=mol
       self.atom=atom 


    def pol_ring(self, len_polymer=2, FF="MMFF",
                 iters=100, shift=1.25): 
      """ Function to create a polymer with ring structure (circular)

      Parameters
      -----------

      mol : rdkit.Chem.rdchem.Mol
            RDKit Mol object
 
      atom : str
            The placeholder atom to combine the molecules and form a new monomer

      len_polymer: int
         Extension of the polymer

      FF: str
         Selected FF to perform a relaxation
 
      iter: int
         Number of iterations to perform a FF geometry optimisation.

      shift: float
         User defined shift for spacing the monomers using the LP function.

      Return
      -------

      pol_ring : rdkit.Chem.rdchem.Mol
           RDKit Mol object
     
      """
      mol=self.mol
      atom=self.atom
   
      new=Lp(mol,str(atom),int(len_polymer), float(shift)).proto_polymer()
      AllChem.EmbedMolecule(new, useRandomCoords=True)

      newMol=self.create_ring(new, str(atom))
      pol_ring=self.check_proto(newMol, str(FF), int(iters))

      return pol_ring
       
    def create_ring(self, mol, atom):
      """Function to create a circular polymer based on a RDKit molecular object.

      Parameters
      ----------

      mol : rdkit.Chem.rdchem.Mol
            RDKit Mol object
 
      atom : str
            The placeholder atom to combine the molecules and form a new monomer
    
      Return
      -------

      final : rdkit.Chem.rdchem.Mol
           RDKit Mol object

      """
      new_bond=seek_plhold(mol, str(atom))

      # New way to flatten a list using sum function
      atoms=sum(new_bond, [])

      rwmol = Chem.RWMol(mol)
      rwmol.AddBond(atoms[0][1], atoms[1][1],
                     rdkit.Chem.rdchem.BondType.SINGLE)

      rwmol.RemoveAtom(atoms[1][0])
      rwmol.RemoveAtom(atoms[0][0])

      final=rwmol.GetMol()
      Chem.SanitizeMol(final)

      return AllChem.AssignBondOrdersFromTemplate(final, final)
   

    def check_proto(self, mol, FF="MMFF", iters=100):
      """Function to check the proto polymer creation.


      Parameters
      -----------
      
      mol: rdkit.Chem.rdchem.Mol
            RDKit Mol object
      
      FF:  str
           Selected FF to perform a relaxation
 
      iter: int
           Number of iterations to perform a FF geometry optimisation.

      Return
      -------

      newMol_H : rdkit.Chem.rdchem.Mol
            RDKit Mol object

      """
      
      newMol_H = Chem.AddHs(mol, addCoords=True)    
      AllChem.EmbedMolecule(mol, useRandomCoords=True)

      # Relaxation functions from utils in linear_polymer
      if FF == "MMFF":
         MMFF_rel(newMol_H, iters)

      else:
         UFF_rel(newMol_H, iters)
    
      return newMol_H


