from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def no_swap(mol, iter_ff, FF="MMFF"):
       """Function to sanitize a molecule with Hydrogens and the user defined atomic place holder.

       Parameters
       -----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

       FF: str
         Selected FF to perform a relaxation
 
       iter_ff: int
         Number of iterations to perform a FF geometry optimization.

       Returns
       --------
       newMol_H : rdkit.Chem.rdchem.Mol
             RDKit Mol object
       
       """
       
       newMol = AllChem.AssignBondOrdersFromTemplate(mol, mol)
       newMol_H = Chem.AddHs(newMol, addCoords=True)

       # This for dealing with big polymers
       Chem.SanitizeMol(newMol_H)
       AllChem.EmbedMolecule(newMol_H, useRandomCoords=True)

       if FF == "MMFF":
           MMFF_rel(newMol_H,iter_ff)

       else:
           UFF_rel(newMol_H, iter_ff)
           

       return newMol_H       


def swap_hyd(mol, iter_ff, atom, FF="MMFF"):
       """Function to swap atomic place holders to Hydrogen atoms.

       Parameters
       -----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

       iter_ff: int
          Number of iterations to perform a FF geometry optimization.

       atom : str
          The placeholder atom to combine the molecules and form a new monomer.

       FF: str
          Selected FF to perform a geometry optimization.
 
       Returns
       ---------

       newMol_H : rdkit.Chem.rdchem.Mol
             RDKit Mol object
       """
       
       for atoms in mol.GetAtoms():
            if atoms.GetSymbol() == str(atom):
                 mol.GetAtomWithIdx(atoms.GetIdx()).SetAtomicNum(1)

       newMol = AllChem.AssignBondOrdersFromTemplate(mol, mol)
       newMol_H = Chem.AddHs(newMol, addCoords=True)

       # This for dealing with big polymers
       Chem.SanitizeMol(newMol_H)
       AllChem.EmbedMolecule(newMol_H, useRandomCoords=True)

       if FF == "MMFF":
           MMFF_rel(newMol_H,iter_ff)

       else:
           UFF_rel(newMol_H, iter_ff)
           

       return newMol_H
       
       
       
def MMFF_rel(mol, iter_ff, vdw_par=0.001):
    """Function to employ a MMFF molecular mechanics FF.
    

    Parameters
    -----------

    mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

    iter_ff: int
          Number of iterations to perform a FF geometry optimization.
          
    vdw_par: float
          Extension of the vdW interaction range. 

    Returns
    ---------

    mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
    """

    AllChem.MMFFOptimizeMolecule(mol, nonBondedThresh=float(vdw_par), maxIters=int(iter_ff))

    return mol

def UFF_rel(mol, iter_ff, vdw_par=0.001):
    """Function to employ an UFF molecular mechanics FF.

    Parameters
    -----------

    mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

    iter_ff: int
          Number of iterations to perform a FF geometry optimization.
          
    vdw_par: float
          Extension of the vdW interaction range. 
    

    Returns
    ---------

    mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
    """
    AllChem.UFFOptimizeMolecule(mol, vdwThresh=float(vdw_par), maxIters=int(iter_ff))

    return mol
