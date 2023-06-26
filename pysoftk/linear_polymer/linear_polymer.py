from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem.rdMolTransforms import *
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


import numpy as np
from pysoftk.tools.utils_func import *
from pysoftk.tools.utils_ob   import *
from pysoftk.tools.utils_rdkit   import *

from openbabel import openbabel as ob
from openbabel import pybel as pb


class Lp:
    """
    A class for creating a linear polymer from given RdKit molecules.

    Examples
    ---------


    Note
    -----

    RDKit package must be installed.
    """
    
    __slots__ = ['mol', 'atom', 'n_copies', 'shift']

    def __init__(self, mol, atom, n_copies, shift):
       """
       Initialize this class.
          
       Parameters
       -----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
 
       atom : str
          The placeholder atom to combine the molecules and form a new monomer

       n_copies: int
          Number of copies to be created for the provided mol object. 

       shift: float
          X-axis shift to translate the super_monomer object.
       """
       
       self.mol = mol
       self.atom = atom
       self.n_copies = n_copies
       self.shift = float(1.25) if shift is None else float(shift)

    def max_dist_mol(self):
        """Returns the maximum distance between atoms from an RDKit Mol Object.

        Return
        -------

        np.float
           maximum value from a Molecule Bound Matrix 
        """
        mol=self.mol
        bm = molDG.GetMoleculeBoundsMatrix(mol)    
        return np.amax(bm)


    def x_shift(self):
        """Re-calibrate distance between monomers
     
        Return
        -------

        shift_final : 'float'
              Value to translates the unit repetitions 
        """
        
        shift=self.shift
        shift_final=float(self.max_dist_mol())-float(shift)
        
        return shift_final

    def _copy_mol(self):
       """Function to replicate super_monomers.

       Parameters
       -----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

       Return
       --------

       fragments : `list`
          List of RDKit Mol objects
       """    

       mol=self.mol
       CanonicalizeConformer(mol.GetConformer())
       
       n_copies=self.n_copies
       
       fragments=[mol for _ in range(int(n_copies))]

       return fragments


    def _polimerisation(self, fragments):
        """Function to produce a polymer in a recursive manner. 

        Parameters
        -----------
       
        fragments: list 
            A list of RDKit objects.

        Return
        --------

        outmol : rdkit.Chem.rdchem.Mol
             RDKit Mol object

        """
        
        x_offset = self.x_shift()
        
        outmol=fragments[0]
        for idx, values in enumerate(fragments[1:]):
            outmol = Chem.CombineMols(outmol,values,
                                      offset=Point3D(x_offset*(idx+1),0.0, 0.0))
            
        order=Chem.CanonicalRankAtoms(outmol, includeChirality=True)
        mol_ordered=Chem.RenumberAtoms(outmol, list(order))

        return outmol

    def _bond_conn(self, outmol):
        """Function to peruse the bonds and connections of place-hold 
        atom within a super_monomer object.
 
        Parameters
        -----------

        outmol : rdkit.Chem.rdchem.Mol
           RDKit Mol object

        Return
        --------

        Tuple : `list`
          A tuple of lists containing connections and 
          neighbours of the place-holder atom.
        
        """
        atom=self.atom

        bonds=atom_neigh(outmol, str(atom))
        conn_bonds=[b for a,b in bonds][1:-1]

        erase_br=[a for a,b in bonds]
        all_conn=list(zip(conn_bonds[::2], conn_bonds[1::2]))

        return all_conn, erase_br
    
    def proto_polymer(self):
       """Function to create a linear polymer from a given 
          super_monomer object.

       Returns
       --------

       newMol_H : rdkit.Chem.rdchem.Mol
            RDKit Mol object

       """
       atom=self.atom
       mol=self.mol
       n_copies=self.n_copies
       
       fragments=self._copy_mol()
       outmol=self._polimerisation(fragments)
       all_conn, erase_br=self._bond_conn(outmol)   

       rwmol = Chem.RWMol(outmol)
       for ini, fin in all_conn:
          rwmol.AddBond(ini, fin, Chem.BondType.SINGLE)

       for i in sorted(erase_br[1:-1], key=None, reverse=True):
          rwmol.RemoveAtom(i)

       mol3=rwmol.GetMol()
       Chem.SanitizeMol(mol3)
       
       mol4=Chem.AddHs(mol3, addCoords=True)
       
       return mol4

    def linear_polymer(self, FF="MMFF94", iter_ff=350, rot_steps=125):
       """Function to create a linear polymer from a 
          given super_monomer object.

       Parameters
       -----------

       FF : str
          Selected Force-Field from RDKit options, i.e, 
          UFF or MMFF.
         
       iter_ff: int
          The maximum number of iterations used for the FF.


       Return
       --------

       newMol_H : rdkit.Chem.rdchem.Mol
            RDKit Mol object

       """

       mol=self.proto_polymer()
       atom=self.atom

       mol1=remove_plcholder(mol, atom)

       #Using PDB object to preserve bond information
       last_rdkit=Chem.MolToPDBBlock(mol1)
       mol_new=pb.readstring('pdb', last_rdkit)

       last_mol=ff_ob_relaxation(mol_new, FF, int(iter_ff))
       rot_mol=rotor_opt(last_mol, FF, int(rot_steps))

       return rot_mol
