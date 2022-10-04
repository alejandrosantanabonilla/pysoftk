from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def swap_hyd(mol, iter_ff, atom, FF="MMFF"):
       """
       Swap Hydrogens
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
       
       
       
def MMFF_rel(mol, iter_ff):
    """
    Function to employ a MMFF molecular mechanics FF.
    
    """

    AllChem.MMFFOptimizeMolecule(mol, maxIters=int(iter_ff))

    return mol

def UFF_rel(mol, iter_ff):
    """
    Function to employ an UFF molecular mechanics FF.
    
    """
    AllChem.UFFOptimizeMolecule(mol, maxIters=int(iter_ff))

    return mol
