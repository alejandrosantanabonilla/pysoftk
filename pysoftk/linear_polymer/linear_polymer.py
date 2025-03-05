import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem.rdMolTransforms import *
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

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

    def copy_mol(self):
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


    def polimerisation(self, fragments):
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

    def bond_conn(self, outmol):
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
       
       fragments=self.copy_mol()
       outmol=self.polimerisation(fragments)
       all_conn, erase_br=self.bond_conn(outmol)   

       rwmol = Chem.RWMol(outmol)
       for ini, fin in all_conn:
          rwmol.AddBond(ini, fin, Chem.BondType.SINGLE)

       for i in sorted(erase_br[1:-1], key=None, reverse=True):
          rwmol.RemoveAtom(i)

       mol3=rwmol.GetMol()
       Chem.SanitizeMol(mol3)
       
       mol4=Chem.AddHs(mol3, addCoords=True)
       
       return mol4

    def linear_polymer(self, force_field="MMFF", relax_iterations=350, rot_steps=125, no_att=True):
        """Function to create a linear polymer from a 
        given super_monomer object.


        Parameters
        -----------

        force_field : str
           Selected Force-Field from RDKit options, i.e., 
           UFF, MMFF, or MMFF94.

        relax_iterations: int
           The maximum number of iterations used for the FF.

        rot_steps: int
           Number of steps for rotor optimization.

        no_att: bool, optional
           If True, remove placeholder atoms. Defaults to True.

        Return
        --------

        newMol_H : rdkit.Chem.rdchem.Mol
           RDKit Mol object

        """

        mol = self.proto_polymer()
        atom = self.atom

        if no_att:
            mol1 = remove_plcholder(mol, atom)
        else:
            mol1 = mol  # Use the original mol if mol1 is not provided

        # Using PDB object to preserve bond information
        last_rdkit = Chem.MolToPDBBlock(mol1)
        mol_new = pb.readstring('pdb', last_rdkit)

        # Validate force field:
        valid_force_fields = ("MMFF", "UFF", "MMFF94")
        if force_field not in valid_force_fields:
            raise ValueError(f"Invalid force field: {force_field}. Valid options are: {valid_force_fields}")

        # Automatically change ff if necessary:
        if force_field == "MMFF":
            force_field = "MMFF94"  # Change to default MMFF94 for MMFF

        # Validate and convert iterations and steps to integers:
        try:
            relax_iterations = int(relax_iterations)
            rot_steps = int(rot_steps)
        except ValueError:
            raise ValueError("relax_iterations and rot_steps must be integers.")

        # Relaxation and optimization:
        last_mol = ff_ob_relaxation(mol_new, force_field, relax_iterations)
        rot_mol = rotor_opt(last_mol, force_field, rot_steps)

        return rot_mol

class Lpr:
    """
    Generates a recursive SMILES string with a stopping condition and a final replacement.
    """
    __slots__ = ['mol', 'replacements', 'max_repetitions', 'final_replacement']

    def __init__(self, mol, replacements, max_repetitions=10, final_replacement="*"):
        """
        Initializes the Lpr object.

        Args:
            smiles: The base SMILES string with a placeholder.
            replacements: A dictionary of replacement strings.
            max_iterations: The maximum number of recursive substitutions.
            final_replacement: The SMILES string or SMART to replace the placeholder after recursion.
        """
        self.mol = mol
        self.replacements = replacements
        self.max_repetitions = max_repetitions
        self.final_replacement = final_replacement

    def generate_recursive_smiles(self, force_field="MMFF", relax_iterations=350, rot_steps=125):
        """
        Performs the recursive SMILES generation.

        Returns:
            A fully resolved SMILES string.
        """
        smiles = self.mol
        for _ in range(self.max_repetitions-1):
            try:
                new_smiles = smiles.format(**self.replacements)
            except KeyError as e:
                print(f"KeyError: {e} in replacements. Please check your replacements dictionary and input SMILES string.")
                return smiles # Return the most recent processed string, or handle as needed.

            if new_smiles == smiles:  # Stop if the SMILES string doesn't change
                break
            smiles = new_smiles

        # Replace the placeholder with the final_replacement
        smiles = smiles.replace("{R}", self.final_replacement)

        mol_new=pb.readstring('smiles', smiles)
        mol_new.make3D()
        
        # Validate force field:
        valid_force_fields = ("MMFF", "UFF", "MMFF94")
        if force_field not in valid_force_fields:
            raise ValueError(f"Invalid force field: {force_field}. Valid options are: {valid_force_fields}")

        # Automatically change ff if necessary:
        if force_field == "MMFF":
            force_field = "MMFF94"  # Change to default MMFF94 for MMFF

        # Validate and convert iterations and steps to integers:
        try:
            relax_iterations = int(relax_iterations)
            rot_steps = int(rot_steps)
        except ValueError:
            raise ValueError("relax_iterations and rot_steps must be integers.")

        # Relaxation and optimization:
        last_mol = ff_ob_relaxation(mol_new, force_field, relax_iterations)
        rot_mol = rotor_opt(last_mol, force_field, rot_steps)
        rel_mol = ff_ob_relaxation(rot_mol, force_field, relax_iterations)
        
        return rel_mol
