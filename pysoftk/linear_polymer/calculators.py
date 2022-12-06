import os
from os import walk
from os.path import join
import sys
import subprocess as sp
import shutil, errno

class Opt(object):
    """ Geometrical optimization tools using GFN-XTB2 and PySCF as engines.


    Note
    -----

    This class requires a GFN-XTB executable and/or a PYSCF semiempirical installation to be used.
    """

    def __init__(self, xyz_file):
        """Set up cartesian coordinates
        
        Parameters
        ----------

        xyz_file : str
            Name of the external cartesian coordinate file used for geometry optimization. 

        """
        self.xyz_file = xyz_file

    def pyscf_semi(self, steps):
        """Function invoking PYSCF semiempirical implementation for geometry optimization.

        Parameters
        ----------

        steps : int, optional
           Number of steps which PYSCF code will perform the optimization routine.


        Return
        -------

        None :
            Print the file pyscf_final.xyz reporting 
            the results from a PYSCF calculation.
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
    """ Print a PYSCF object into cartesian 
        coordinates.
    """

    def xyz(self, mol):
        """Function to print a file with Cartesian coordinates from a PYSCF geometry optimization calculation.

        Parameters
        ----------
       
        mol : pyscf.gto.mole.Mole
            A PySCF gto mole Mole object.

        Return
        ------

        pysfc_final: str
          An xyz file containing the optimized coordinates of the provided PYSCF Mole object.
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



