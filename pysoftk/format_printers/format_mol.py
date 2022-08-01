from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG

class Fmt(object):
    """Class printing RDkit Mol Object in 
       different formats
          
    Returns
    -------
    file : str
            A file with the provided molecule
            in a chosen format.

    Note
    ----
    This class requires MDanalysis to be installed. 
    """
    
    def __init__(self, mol):
       """ Initialise the class Fmt

       Parameters
       ----------
       mol : Chem.rdchem.Mol
           an RDKit Mol object Molecule to be 
           printed in a given format.
       """
       self.mol = mol
    
    def xyz_print(self, output_name):
        """Function to print in XYZ format
        
        Parameters
        ----------
        output_name : str 
           A provided name for an output file.
     
        Raises
        ------
        ValueError:
            If not an RDKit Mol Object raises 
            TypeError as invalid molecule.
        """
        mol=self.mol
        
        try:
          isinstance(mol,Chem.rdchem.Mol)
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))

        else:
          Chem.MolToXYZFile(mol, str(output_name))

    def pdb_print(self, output_name):
        """Function to print in PDB format

        Parameters
        ----------
        output_name : str 
           A provided name for a file.
 
        Raises
        ------
        ValueError:
            If not an RDKit Mol Object raises 
            TypeError as invalid molecule.
        """
        mol=self.mol
        try:
          isinstance(mol,Chem.rdchem.Mol)
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))
        else:
          Chem.MolToPDBFile(mol,str(output_name))
            
          
    def mol_print(self, output_name):
        """Function to print in MOL format

        Parameters
        ----------
        output_name : str 
           A provided name for a file.
 
        Raises
        ------
        ValueError:
            If not an RDKit Mol Object raises 
            TypeError as invalid molecule.
        """
        mol=self.mol
        try:
          isinstance(mol,Chem.rdchem.Mol)
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))
        else:
          Chem.MolToMolFile(mol,str(output_name))


