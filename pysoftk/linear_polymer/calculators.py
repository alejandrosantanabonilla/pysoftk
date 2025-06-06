import os
from os import walk
from os.path import join
import sys
import subprocess as sp
import shutil, errno

class Opt(object):
    """
    Provides tools for performing **geometrical optimizations** of molecular structures.
    This class leverages external computational chemistry engines, specifically
    **GFN-XTB2** and **PySCF** (with its semiempirical implementations), to minimize
    the energy of a given molecular geometry.

    The primary use case for this class is to find stable 3D conformations of molecules
    by iteratively adjusting atomic positions until a local energy minimum is reached.

    Note
    -----
    Successful execution of the optimization methods within this class requires
    the prior installation and correct configuration of the chosen computational
    engine:
    * **GFN-XTB**: The `xtb` executable must be accessible in your system's PATH.
        You can typically install it via conda or by compiling from source.
    * **PySCF Semiempirical**: The `pyscf` Python package, including its
        semiempirical modules, must be installed in your Python environment.
        Installation is usually done via pip (`pip install pyscf`).
    """

    def __init__(self, xyz_file):
        """
        Initializes the `Opt` class by setting the path to the input
        Cartesian coordinate file. This file serves as the starting
        geometry for all subsequent optimization routines.

        Parameters
        ----------
        xyz_file : str
            The **name or path** of the external file containing the
            Cartesian coordinates of the molecular structure to be optimized.
            This file should typically be in the `.xyz` format.
        """
        
        self.xyz_file = xyz_file

    def pyscf_semi(self, steps):
        """
        Performs a **geometry optimization** using the **PySCF semiempirical**
        implementation, specifically the **MINDO/3** method.

        This function sets up a PySCF molecular object from the provided XYZ file,
        configures a semiempirical calculation, and then applies a geometry
        optimization algorithm (Berny solver) to find the equilibrium structure.
        The optimized coordinates are then printed to an XYZ file.

        Parameters
        ----------
        steps : int, optional
            The **maximum number of steps** the PySCF optimization routine
            will attempt to converge the geometry. The optimization may
            terminate earlier if convergence criteria are met.

        Returns
        -------
        None
            This method does not return any value. Instead, it **generates an output file**
            named `pyscf_final.xyz` in the current working directory, which
            contains the optimized Cartesian coordinates of the molecule.
            A confirmation message will be printed to the console upon successful
            file creation.
        """
        
        import pyscf
        from pyscf import gto
        from pyscf import semiempirical
        from pyscf.geomopt.berny_solver import optimize

        mol=gto.M(atom=str(self.xyz_file))
        mf=semiempirical.RMINDO3(mol)
        mol_eq=optimize(mf,maxsteps=int(steps))
        Pyscf_print().xyz(mol_eq)

class Pyscf_print(object):
    """
    A utility class designed to **convert and save** the results of a PySCF
    geometry optimization into a standard **Cartesian coordinate (XYZ)** file format.

    This class acts as a bridge between PySCF's internal molecular object
    representation and a widely recognized format for molecular structures,
    making the optimized geometries easily shareable and viewable with
    molecular visualization software.
    """

    def xyz(self, mol):
        """
        Writes the **optimized Cartesian coordinates** from a PySCF molecular
        object to a new `.xyz` file named `pyscf_final.xyz`.

        This function takes the atomic symbols and their optimized coordinates
        from a PySCF `Mole` object, converts the coordinates from Bohr to Angstroms,
        and formats them into the XYZ file standard.

        Parameters
        ----------
        mol : pyscf.gto.mole.Mole
            A **PySCF `gto.Mole` object** containing the optimized molecular geometry,
            typically obtained from a PySCF geometry optimization calculation.
            This object holds information about atomic symbols and their Cartesian coordinates.

        Returns
        -------
        None
            This method does not return any value. It **creates a file** named
            `pyscf_final.xyz` in the current working directory, which contains
            the optimized molecular structure. A confirmation message
            `"File pyscf_final.xyz has been created."` is printed to standard output.
        """
        
        bohr2angs=0.529177
        coords = bohr2angs*mol.atom_coords()
        with open("pyscf_final.xyz", "w") as myfile:
            myfile.write('{}'.format(int(len(coords))))
            myfile.write('\n')
            myfile.write('\n')
            for i in range(len(mol.atom)):
                myfile.write('{} {:.8f} {:.8f} {:.8f}\n'.format(str(mol.atom[i][0]),
                                                                 float(coords[i][0]),
                                                                 float(coords[i][1]),
                                                                 float(coords[i][2])))
            myfile.write('\n')
        print ("File pyscf_final.xyz has been created.")
