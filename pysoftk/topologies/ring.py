import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import copy

from pysoftk.topologies.diblock import *
from pysoftk.linear_polymer.linear_polymer import *

from pysoftk.tools.utils_func import *
from pysoftk.tools.utils_rdkit import *

class Rn:
    """A class for creating a circular (ring-shaped) polymer 
       from given RDKit molecules.

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
                 iter_ff=100, shift=1.25): 
      """ Function to create a polymer with ring structure (circular)

      Parameters
      -----------

      mol : rdkit.Chem.rdchem.Mol
            RDKit Mol object
 
      atom : str
            The placeholder atom to combine the molecules and 
            form a new monomer

      len_polymer: int
         Extension of the polymer

      FF: str
         Selected FF to perform a relaxation
 
      iter_ff: int
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
   
      #new=Lp(mol,str(atom),int(len_polymer), float(shift)).proto_polymer()
      patt=Lp(mol,str(atom),int(len_polymer),float(shift)).proto_polymer()
      
      newMol=self.create_ring(patt, str(atom))
      pol_ring=self.check_proto(newMol, str(FF), int(iter_ff))

      return pol_ring
       
    def create_ring(self, mol, atom):
      """Function to create a circular polymer based on a RDKit 
         molecular object.

      Parameters
      ----------

      mol : rdkit.Chem.rdchem.Mol
            RDKit Mol object
 
      atom : str
            The placeholder atom to combine the molecules and 
            form a new monomer
    
      Return
      -------

      final : rdkit.Chem.rdchem.Mol
           RDKit Mol object

      """
     
      atoms=atom_neigh(mol, str(atom))
    
      rwmol = Chem.RWMol(mol)
      rwmol.AddBond(atoms[0][1], atoms[1][1],
                     rdkit.Chem.rdchem.BondType.SINGLE)

      rwmol.RemoveAtom(atoms[1][0])
      rwmol.RemoveAtom(atoms[0][0])

      final=rwmol.GetMol()
      Chem.SanitizeMol(final)

      newMol_H = Chem.AddHs(final, addCoords=True)  

      return AllChem.AssignBondOrdersFromTemplate(newMol_H, newMol_H)
   

    def check_proto(self, mol, FF="MMFF", iter_ff=100):
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
      last_rdkit=Chem.MolToPDBBlock(mol)

      mol_new=pb.readstring('pdb', last_rdkit)

      # Relaxation functions from utils in linear_polymer
      opt_mol=ff_ob_relaxation(mol_new, iter_ff=100, ff_thr=1.0e-6)


      return opt_mol

      



