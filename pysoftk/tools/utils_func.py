import os

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def get_file_extension(path):
    """Gets the extension of a file.

    Args:
        path (str): The path to the file.

    Returns:
        str: The extension of the file.
    """

    file_name, file_extension = os.path.splitext(path)
    
    return file_name,file_extension


def pattern_recon(pattern):
    """Function to find unique values from a string pattern.
    
    Parameters
    ------------  

    pattern: str
        Variable containing the pattern written in a string.
    
    Returns
    ---------

    result: list    
        Returns the unique values of the pattern in alphabetic order.
    """
    return sorted(list(set(pattern)), key=lambda c:c.upper())

def pattern_repl(pattern, tup_repl):
   """Function to perform a pattern replacement using tuples generated 
      by pattern and molecular names.
   
   Parameters
   ------------  

   pattern: str
        Variable containing the pattern written in a string

   
   Return
   ---------

   return: list    
        Returns the names of the molecules as provided 
        by the string pattern.
   """
   
   for r in tup_repl:
      pattern=pattern.replace(*r)

   valid_seq=pattern.split("+")
   return list(filter(None, valid_seq))


def pattern_mol_seq(mols, pattern):
    """Function to create a list of molecules based on a user 
       provided pattern.

    Parameters
    ------------
  
    pattern: class.str
         A string value containing the desired pattern.

    mols: class.list
        List containing the names of the molecules to form a 
        polymer with a given pattern.

    Returns
    ---------

    return: class.list
      List with the correct sequence of molecular fragments 
      to be merged. 
  
    """
    unq_val=pattern_recon(pattern)
    tup_repl=tuple(zip(unq_val,mols))

    return pattern_repl(pattern, tup_repl)

def count_plholder(mol, atom):
    """Function to count the number of user defined 
       atomic placeholders.

    Parameters
    -----------

    mol: rdkit.Chem.rdchem.Mol
          RDKit Mol object

    atom: class.str
         The placeholder atom
    
    Returns
    -------

    num_br: class.int
        Number of atomic place holders inside a molecule.

    """
    num_br=0
    
    for atoms in mol.GetAtoms():
         if atoms.GetSymbol() == str(atom):
                num_br=num_br+1

    return num_br

def atom_neigh(mol, atom):
   """Find the neighbours from an atoms inside a molecule. 

   Parameters
   -----------

   mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

   atom : str
         The placeholder atom

   Returns
   --------

   neigh: list
         List of neighbors

   """

   neigh=[]
   rwmol = Chem.RWMol(mol)
   for atoms in rwmol.GetAtoms():
       if atoms.GetSymbol() == str(atom):
             neigh.append([(atoms.GetIdx(),nbr.GetIdx())
                           for nbr in atoms.GetNeighbors()])

   return sum(neigh, [])

def tuple_bonds(lst_atm_neigh):
   """Function to seek atomic neighbours

   Parameters
   -----------

   lst_atm_neigh: class.list
                  List containing atoms for getting neighbors.

   Return
   -------

   return: class.list
           A list of tuples with the neighbors of the parsed atoms. 
   
   """

   bond_idx=[]
   for i in range(1,len(lst_atm_neigh)-1):
       a,b = lst_atm_neigh[i]
       bond_idx.append(b)

   return [(bond_idx[i],bond_idx[i+1])
           for i in range(0,len(bond_idx),2)]

def create_pol(mol, atom, tpb):
    """Function to create polymers.

    Parameters
    -----------

    mol: rdkit.Chem.rdchem.Mol
          RDKit Mol object

    atom: str
          The placeholder atom

    
    tpb: list
         List of tuples with atoms to be connected


    Returns
    --------

    rmwol1: rdkit.Chem.rdchem.Mol
          RDKit Mol object
    

    """

    rwmol1 = Chem.RWMol(mol)
    for ini, fin in tpb:
        rwmol1.AddBond(ini, fin, Chem.BondType.SINGLE)

    #Reading Place holder as SMILES
    plc=Chem.MolFromSmiles(str(atom))
    new_conn=sum(rwmol1.GetSubstructMatches(plc),())

    for i in sorted(new_conn[1:-1], key=None, reverse=True):
         rwmol1.RemoveAtom(i)
         
    return rwmol1
