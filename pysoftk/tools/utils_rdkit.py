from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from openbabel import openbabel as ob
from openbabel import pybel as pb

def no_swap(mol, iter_ff, force_field="MMFF"):
       """Function to sanitize a molecule with Hydrogens 
          and the user defined atomic place holder.

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

       if force_field == "MMFF":
           MMFF_rel(newMol_H,iter_ff)

       else:
           UFF_rel(newMol_H, iter_ff)
           

       return newMol_H       


def swap_hyd(mol, iter_ff, atom, force_field="MMFF"):
       """Function to swap atomic place holders to Hydrogen atoms.

       Parameters
       -----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

       iter_ff: int
          Number of iterations to perform a FF geometry optimization.

       atom : str
          The placeholder atom to combine the molecules and form a new monomer.

       force_field: str
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

       if force_field == "MMFF":
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


def plc_holder(mol, atom):
    """Function Seeking for a specific placeholder atom.

    Parameters
    -----------

    mol: rdkit.Chem.rdchem.Mol
         RDKit Mol object

    atom : str
         The placeholder atom


    Returns
    --------

    return: list
         List of neighbors from the place holder atom.
    
    """
    
    new_bond=[]

    for atoms in mol.GetAtoms():
        if atoms.GetSymbol() == str(atom):
            new_bond.append([[ atoms.GetIdx(), nbr.GetIdx()]
                             for nbr in atoms.GetNeighbors()])


    return new_bond

def remove_plcholder(mol, atom):
    """Function that seeks for a place holder atom 
       and replace it with a Hydrogen atom.

      Parameters
      ----------
       
      mol: rdkit.Chem.Mol
          RDKit Mol object


     Return
     --------

     None:
          RDKit Mol object
    

    """
    for atoms in mol.GetAtoms():
        if atoms.GetSymbol() == str(atom):
              mol.GetAtomWithIdx(atoms.GetIdx()).SetAtomicNum(1)

    Chem.SanitizeMol(mol)

    return mol

def etkdgv3_energies(mol, num_conf=1):
   """Calculate molecular configurations using the RDKit-ETKDG3 method.
      
      Parameters
      ----------
       
      mol: rdkit.Chem.Mol
          RDKit Mol object

      num_conf: int
          The number of configurations requested to be computed.

     
      Return
      --------

      datapoint: rdkit.Chem.rdistGeom.EmbedMultipleConfs
            RDKit Mol object

   """
   
   AllChem.EmbedMolecule(mol, useRandomCoords=True)
   AllChem.MMFFOptimizeMolecule(mol)
     
   m=Chem.Mol(mol)
      
   ps=AllChem.ETKDGv3()
   ps.randomSeed=0xf00d
   ps.RandomCoords=False
     
   cids=AllChem.EmbedMultipleConfs(m,int(num_conf),ps)
   mp=AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
   AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0,
                                     mmffVariant='MMFF94s')

   energies=[]
   for cid in cids:
      ff = AllChem.MMFFGetMoleculeForceField(m, mp, confId=cid)
      energies.append(ff.CalcEnergy())

   return m, cids, energies 
