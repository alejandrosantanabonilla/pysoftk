from openbabel import openbabel as ob
from openbabel import pybel as pb

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG

from pysoftk.tools.utils_func import *

class Fmt(object):
    """Class printing RDkit Mol Object in different formats
          
    Returns
    -------

    file : str
            A file with the provided molecule in a chosen format.

    """
    
    def __init__(self, mol):
       """ Initialise the class Fmt

       Parameters
       ----------

       mol : Chem.rdchem.Mol
           An RDKit Mol object Molecule to be printed in a given format.
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
            If not an RDKit Mol Object raises TypeError as invalid molecule.
        """
        mol=self.mol
        
        try:
           if isinstance(mol,Chem.rdchem.Mol):
               Chem.MolToXYZFile(mol, str(output_name))

           else:
               mol.OBMol.PerceiveBondOrders()
               mol.OBMol.AddHydrogens()
               mol.write('xyz',str(output_name), overwrite=True)
               
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))

        
          

    def pdb_print(self, output_name):
        """Function to print in PDB format

        Parameters
        ----------

        output_name : str 
           A provided name for a file.
 
        Raises
        ------

        ValueError:
            If not an RDKit Mol Object raises TypeError as invalid molecule.
        """
        mol=self.mol
        
        try:
           if isinstance(mol,Chem.rdchem.Mol):
               Chem.MolToPDBFile(mol,str(output_name))

           else:
               mol.OBMol.PerceiveBondOrders()
               mol.OBMol.AddHydrogens()
               mol.write('pdb',str(output_name), overwrite=True)
               
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))
            
          
    def mol_print(self, output_name):
        """Function to print in MOL format

        Parameters
        ----------

        output_name : str 
           A provided name for a file.
 
        Raises
        ------
        ValueError:
            If not an RDKit Mol Object raises TypeError as invalid molecule.
        """
        mol=self.mol

        try:
           if isinstance(mol,Chem.rdchem.Mol):
               Chem.MolToMolFile(mol,str(output_name))

           else:
               mol.OBMol.PerceiveBondOrders()
               mol.OBMol.AddHydrogens()
               mol.write('mol',str(output_name), overwrite=True)
               
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))


    def format_print(self, output_name):
        """Function to print a pysoftk.object into an user 
        provided output_format.

        Parameters
        ----------

        fmt: class.str
           An user-provided valid format in which the pysoftk.object
           can be printed. Default is MOL2 format.
 
        output_name : class.str 
           A provided name for a file. Default is molecule.mol2 
 
        """

        mol=self.mol
        ext=output_name.split('.')[-1]
        
        try:
           if isinstance(mol,Chem.rdchem.Mol):
              last_rdkit=Chem.MolToPDBBlock(mol)
              mol2=pb.readstring('pdb', last_rdkit)
              mol2.OBMol.PerceiveBondOrders()
              mol2.OBMol.AddHydrogens()
              mol2.write(str(ext),str(output_name), overwrite=True)

           else:
               mol.OBMol.PerceiveBondOrders()
               mol.OBMol.AddHydrogens()
               mol.write(str(ext),str(output_name), overwrite=True)
               
        except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))
     

class Cnv(object):
    """Class for converting files from X to Y format. 
          
    Returns
    -------

    file : str
            A file with the provided molecule in a chosen format.
    """

    def __init__(self, file_mol):
       """ Initialise the class Fmt

       Parameters
       ----------

       file_mol : class.str
          User-provided file name. 
       
       """

       self.file_mol = file_mol
       
    def file_in_out(self, fmt="mol2"):
        """Function to convert an user-provided file 
           with X format to an user-provided Y format.

           Parameters
           ------------

           fmt: class.str
             Format requested to convert.

           Return
           -------

           None:
             New file with the provided format.

        """
        file_mol=self.file_mol 
        name, extension=get_file_extension(str(file_mol))
        
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats(str(extension), str(fmt))

        mol = ob.OBMol()
        obConversion.ReadFile(mol, str(file_mol)) 
        obConversion.WriteFile(mol, '.'.join([name,str(fmt)]))
