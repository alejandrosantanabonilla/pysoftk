from openbabel import openbabel as ob
from openbabel import pybel as pb

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG

from pysoftk.tools.utils_func import *

class Fmt(object):
    """
    A class designed for printing RDKit Mol objects into various molecular file formats.

    This class provides methods to convert an RDKit Mol object into common chemical file
    formats such as XYZ, PDB, and MOL, as well as a more general format conversion.

    Returns
    -------
    file : str
        A file containing the provided molecule in the chosen format.

    """
    
    def __init__(self, mol):
        """
        Initializes the Fmt class with an RDKit Mol object.

        Parameters
        ----------
        mol : Chem.rdchem.Mol or pybel.Molecule
            The molecule to be processed. This can be an RDKit Mol object or a
            Pybel Molecule object.
        """
        
        self.mol = mol
    
    def xyz_print(self, output_name):
        """
        Prints the molecule to an XYZ file format.

        This method handles both RDKit Mol objects and Pybel Molecule objects.
        For RDKit Mol objects, it uses `Chem.MolToXYZFile`. For Pybel Molecule
        objects, it ensures bond orders are perceived and hydrogens are added
        before writing the file.

        Parameters
        ----------
        output_name : str
            The desired name for the output XYZ file.

        Raises
        ------
        TypeError
            If the input `mol` is neither an RDKit Mol object nor a Pybel Molecule
            object, indicating an invalid molecule type.
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
            raise TypeError(f"{mol} is an invalid molecule")

    def pdb_print(self, output_name):
        """
        Prints the molecule to a PDB file format.

        This method supports both RDKit Mol objects and Pybel Molecule objects.
        For RDKit Mol objects, it uses `Chem.MolToPDBFile`. For Pybel Molecule
        objects, it ensures bond orders are perceived and hydrogens are added
        before writing the file.

        Parameters
        ----------
        output_name : str
            The desired name for the output PDB file.

        Raises
        ------
        TypeError
            If the input `mol` is neither an RDKit Mol object nor a Pybel Molecule
            object, indicating an invalid molecule type.
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
            raise TypeError(f"{mol} is an invalid molecule")
            
    def mol_print(self, output_name):
        """
        Prints the molecule to a MOL file format.

        This method handles both RDKit Mol objects and Pybel Molecule objects.
        For RDKit Mol objects, it utilizes `Chem.MolToMolFile`. For Pybel Molecule
        objects, it ensures bond orders are perceived and hydrogens are added
        before writing the file.

        Parameters
        ----------
        output_name : str
            The desired name for the output MOL file.

        Raises
        ------
        TypeError
            If the input `mol` is neither an RDKit Mol object nor a Pybel Molecule
            object, indicating an invalid molecule type.
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
            raise TypeError(f"{mol} is an invalid molecule")

    def format_print(self, output_name):
        """
        Prints the molecule into a user-provided file format.

        This function attempts to convert and write the molecule to a file
        based on the extension provided in `output_name`. It handles both RDKit
        Mol objects and Pybel Molecule objects. If an RDKit Mol is provided,
        it's first converted to a PDB block and then read by Pybel to facilitate
        writing to the desired format.

        Parameters
        ----------
        output_name : str
            The desired name for the output file, including the file extension
            (e.g., "molecule.mol2", "compound.smi"). The extension determines
            the output format.

        Raises
        ------
        TypeError
            If the input `mol` is neither an RDKit Mol object nor a Pybel Molecule
            object, indicating an invalid molecule type.
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
            raise TypeError(f"{mol} is an invalid molecule")
    

class Cnv(object):
    """
    A class for converting molecular files from one format to another using OpenBabel.

    This class provides a utility to convert a molecule stored in a file of a
    given format (e.g., PDB, SDF, MOL2) to a different specified format.

    Returns
    -------
    file : str
        A new file containing the provided molecule in the chosen output format.
    """

    def __init__(self, file_mol):
        """
        Initializes the Cnv class with the path to a molecular file.

        Parameters
        ----------
        file_mol : str
            The path to the input molecular file that needs to be converted.
        """

        self.file_mol = file_mol
        
    def file_in_out(self, fmt="mol2"):
        """
        Converts a molecular file from its current format to a specified output format.

        This function leverages OpenBabel to perform the file conversion. It
        automatically detects the input file format from its extension and
        converts it to the format specified by the `fmt` parameter.

        Parameters
        ----------
        fmt : str, optional
            The desired output format for the converted file (e.g., "pdb", "xyz", "sdf").
            Defaults to "mol2".

        Returns
        -------
        None
            This function does not return any value. It creates a new file
            in the same directory as the input file, with the converted format
            and the original filename (before its extension) appended with the
            new format's extension.
        """
        
        file_mol=self.file_mol 
        name, extension=get_file_extension(str(file_mol))
        
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats(str(extension), str(fmt))

        mol = ob.OBMol()
        obConversion.ReadFile(mol, str(file_mol)) 
        obConversion.WriteFile(mol, '.'.join([name,str(fmt)]))
