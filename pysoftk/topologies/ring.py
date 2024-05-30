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

from openbabel import openbabel as ob
from openbabel import pybel as pb

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
<<<<<<< HEAD
     
=======

    def check_proto(self, mol, force_field="MMFF", iter_ff=100):
      """Function to check the proto polymer creation.


      Parameters
      -----------
      
      mol: rdkit.Chem.rdchem.Mol
            RDKit Mol object
      
      force_field:  str
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

      if force_field == "MMFF":
         ff = pb._forcefields["mmff94"]
         # Relaxation functions from utils in linear_polymer
         opt_mol=ff_ob_relaxation(mol_new, iter_ff=int(iter_ff), ff_thr=1.0e-6)

      else:
         ff = pb._forcefields["uff"]
         # Relaxation functions from utils in linear_polymer
         opt_mol=ff_ob_relaxation(mol_new, iter_ff=int(iter_ff), ff_thr=1.0e-6)
          
      return opt_mol
       

>>>>>>> 88f2163339a8017b72d141956fc54ab10ba533f5
    def pol_ring(self, len_polymer=2, force_field="MMFF",
                 relax_iterations=100, shift=1.25, more_iter=10): 
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

      force_field: str
         Selected FF to perform a relaxation
 
      relax_iterations: int
         Number of iterations to perform a FF geometry optimisation.

      shift: float
         User defined shift for spacing the monomers using the LP function.

      more_iter: int
         User defined variable to perform more iterations using the selected FF.

      Return
      -------

      pol_ring : rdkit.Chem.rdchem.Mol
           RDKit Mol object
     
      """
      mol=self.mol
      atom=self.atom


      if force_field not in ("MMFF", "UFF"):
            raise ValueError(f"Invalid force field: {force_field}")
      
      patt=Lp(mol,str(atom),int(len_polymer),float(shift)).linear_polymer(relax_iterations=int(relax_iterations),
                                                                          no_att=False)
      proto_ring=patt.write("pdb")

      ring_rdkit=Chem.MolFromPDBBlock(proto_ring)   
      atoms=atom_neigh(ring_rdkit, str(atom))
      
      rwmol=Chem.RWMol(ring_rdkit)

      rwmol.AddBond(atoms[0][1], atoms[1][1],
                    Chem.BondType.SINGLE)
      
      rwmol.RemoveAtom(atoms[1][0])
      rwmol.RemoveAtom(atoms[0][0])

      final=rwmol.GetMol()
      Chem.SanitizeMol(final)

      newMol_H = Chem.AddHs(final, addCoords=True)  
      AllChem.AssignBondOrdersFromTemplate(newMol_H, newMol_H)

      if force_field == "MMFF":
         apply_force_field(newMol_H, "MMFF", relax_iterations*int(more_iter))

      else:
         apply_force_field(newMol_H, "UFF", relax_iterations*int(more_iter)) 

      pol_ring=check_proto(newMol_H, str(force_field), int(relax_iterations))
         
      return pol_ring
      
     
<<<<<<< HEAD
def apply_force_field(molecule: Chem.Mol, force_field: str, relax_iterations: int) -> None:
    """Applies the specified force field to the molecule."""

    if force_field == "MMFF":
        MMFF_rel(molecule, relax_iterations)

    else:
        UFF_rel(molecule, relax_iterations)
=======
def apply_force_field(molecule: Chem.Mol, force_field: str, iter_ff: int) -> None:
  """Applies the specified force field to the molecule.

  This function takes a RDKit molecule object (`molecule`), a string specifying the force field to use (`force_field`),
  and an integer representing the number of iterations (`iter_ff`). It then applies the chosen force field to the molecule
  using either the MMFF or UFF approach, depending on the `force_field` value.

  Args:
    molecule: A RDKit molecule object representing the molecule for force field application.
    force_field: A string specifying the force field to use, either "MMFF" or "UFF".
    iter_ff: An integer representing the number of iterations to perform during force field minimization.

  Raises:
    ValueError: If an invalid force field name is provided.

  Returns:
    None. The function modifies the input molecule object `molecule` in-place.

  """
  if force_field == "MMFF":
     MMFF_rel(molecule, iter_ff)

  else:
     UFF_rel(molecule, iter_ff)
>>>>>>> 88f2163339a8017b72d141956fc54ab10ba533f5


def check_proto(mol, force_field="MMFF", relax_iterations=100, rot_steps=1, ff_thr=1.0e-6):
   """Function to check the proto polymer creation.


   Parameters
   -----------
      
   mol: rdkit.Chem.rdchem.Mol
        RDKit Mol object
      
   force_field:  str
        Selected FF to perform a relaxation
 
   relax_iterations: int
        Number of iterations to perform a FF geometry optimisation.

   rot_steps: int 
        Number of possible confomers to search after optimisation.

   ff_thr: force field convergence threshold
        Force convergence criteria.

   Return
   -------

   newMol_H : rdkit.Chem.rdchem.Mol
        RDKit Mol object

   """
   last_rdkit=Chem.MolToPDBBlock(mol)

   mol_new=pb.readstring('pdb', last_rdkit)

   if force_field == "MMFF":
      ff = pb._forcefields["mmff94"]
      # Relaxation functions from utils in linear_polymer
      #opt_mol=ff_ob_relaxation(mol_new, relax_iterations=int(relax_iterations), ff_thr=1.0e-6)
      opt_mol=global_opt(mol_new, relax_iterations=int(relax_iterations), rot_steps=int(rot_steps),ff_thr=float(ff_thr))

   else:
      ff = pb._forcefields["uff"]
      # Relaxation functions from utils in linear_polymer
      #opt_mol=ff_ob_relaxation(mol_new, relax_iterations=int(relax_iterations), ff_thr=1.0e-6)
      opt_mol=global_opt(mol_new, relax_iterations=int(relax_iterations), rot_steps=int(rot_steps),ff_thr=float(ff_thr))
          
   return opt_mol
