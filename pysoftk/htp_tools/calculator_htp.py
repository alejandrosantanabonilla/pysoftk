import glob
import subprocess as sp
from   subprocess import call
from concurrent.futures import ThreadPoolExecutor, as_completed
import os, fnmatch
from tqdm import tqdm
import time

class Htp(object):
    """Class enabling High-throughput calculations 
       in different created folders for provided
       molecules

    Attributes 
    ----------
    ending : str
       Extension of the molecular format used to perform a
       search within many directories.
    xtb_path : str
          Path where the xtb executable is located. 
    xyz_geom : str
          Name of the supplied molecule in .xyz format
    init_dir : str
          Path where the calculation is starting (optional)
    num_cores : int, optional
          Number of cores which the GFN-XTB2 code
          will use.
    threshold: str, optional
          Level of geometry optimization required 
          as defined in GFN-XTB2 code. 
        
    Note
    -----
    This class requires tqdm, GFN-XTB and/or PYSCF semiempirical 
    installation to be used. 
    """

    def __init__(self, ending):
      """Parameters to initialize the HTP process

        Parameters
        ----------
        ending : str
           Name of the extension to be browsed in a
           provided folder. 
      """
      self.ending = ending


    @staticmethod   
    def xtb_gfn(xtb_path, xyz_geom, init_dir, num_cores=None, threshold=None):
       """Function invoking GFN-XTB Force-Field implementation for 
           geometry optimization.

       Parameters
       -----------
       xtb_path : str
           Path where the GFN-XTB2 executable is located

       num_cores : int, optional
             Number of cores which the GFN-XTB2 code will use.

       threshold : str, optional
             Level of geometry optimization required 
             as defined in GFN-XTB2 code.
 
       xyz_geom : str 
             Provided geometry to be relaxed.

       init_dir : str
             Root path directory where folders containing 
             the structures are found.
 
       Returns
       -------
       output : str
          File name output.log containing the optimization
          process carried out by GFN-XTB2.
       """
    
       num_cores = int(1) if num_cores is None else int(num_cores)
       threshold = str("crude") if threshold is None else str(threshold)
    
       out_file = os.path.join(init_dir,"output.log")

       cmd = (
              f'{xtb_path} {xyz_geom} '            
              f'--parallel {num_cores} '
              f'--opt {threshold} '
             )

          # Note that sp.call will hold the program until completion
          # of the calculation.

       try:
          with open(out_file, 'w') as result:
            sp.call(
               cmd,
               stdin=sp.PIPE,
               stdout=result,
               stderr=sp.PIPE,
               shell=True,
               cwd=init_dir
              )
       finally:
           pass

        
    @staticmethod   
    def xtb_ff(xtb_path, xyz_geom, init_dir, num_cores=None, threshold=None):
       """Function invoking GFN-XTB Force-Field implementation for geometry optimization.

          Parameters
          ------------
          xtb_path : str
             Path where the GFN-XTB2 executable is 
             located
          num_cores : int, optional
             Number of cores which the GFN-XTB2 code
             will use.
          threshold : str, optional
             Level of geometry optimization required 
             as defined in GFN-XTB2 code. 
          xyz_geom : str 
             Provided geometry to be relaxed.
          init_dir : str
             Root path directory where folders containing 
             the structures are found.
 
         Returns
         -------
          output : str
             File name output.log containing the optimization
             process carried out by GFN-XTB2.
       """
    
       num_cores = int(1) if num_cores is None else int(num_cores)
       threshold = str("crude") if threshold is None else str(threshold)
    
       out_file = os.path.join(init_dir,"output.log")

       cmd = (
              f'{xtb_path} --gfnff {xyz_geom} '            
              f'--parallel {num_cores} '
              f'--opt {threshold} '
             )

          # Note that sp.call will hold the program until completion
          # of the calculation.

       try:
          with open(out_file, 'w') as result:
            sp.call(
               cmd,
               stdin=sp.PIPE,
               stdout=result,
               stderr=sp.PIPE,
               shell=True,
               cwd=init_dir
              )
       finally:
           pass

    @staticmethod
    def pyscf_semi(steps, xyz_geom, save_path):
        """Function invoking PYSCF semiempirical implementation
           for geometry optimization.

        Parameters
        ----------
        steps : int, optional
           Number of steps which PYSCF code
           will perform the optimization routine.

        Results
        -------
        None :
            Print the file pyscf_final.xyz reporting 
            the results from a PYSCF calculation.
        """
        import pyscf
        from pyscf import gto
        from pyscf import semiempirical
        from pyscf.geomopt.berny_solver import optimize
    
        mol=gto.M(atom=str(xyz_geom))
        mf=semiempirical.RMINDO3(mol)
        mol_eq=optimize(mf,maxsteps=int(steps))
        Pyscf_print().xyz(mol_eq, save_path)

        
    def _find_files(self, directory, pattern):
        """Function to find files within a given directory
           based on a provided pattern.

        Parameters
        ----------
        directory : str
           Path of the main directory where the search will
           be carried out. 

        pattern : str
           Provided Alpha-pattern (extension of the molecular
           file) which will be used to perform the search.

        Results
        -------
        None :
            The results of the search as a string

        """
        for root, dirs, files in os.walk(directory):
            for basename in files:
                if fnmatch.fnmatch(basename, pattern):
                    filename = os.path.join(root, basename)
                    yield filename

    # High-throughput Part of the Class

    def _seek_files(self, file_ending):
        """Function to seek files within the current working 
           directory based on a provided pattern.

        Parameters
        ----------
        file_ending : str
           Extension of the provided pattern to perform 
           the search.

        Results
        -------
        gmts : list
           A list of the resulting files within the 
           current working directory with the provided
           extension (ending). 

        """
        seek_file = "".join(("*",".", str(file_ending)))
        gmts = [filename for filename
                in self._find_files(os.getcwd(), str(seek_file))]
        return gmts


    def htp_xtb_ff(self, path_xtb, max_work, cores, threshold="crude"):
       """Function to perform a parallel concurrent calculation
          at the gfn-ff level of theory.

       Parameters
       ----------
       path_xtb : str
          Path where the xtb executable is located. 
       max_work : int
          Number of maximum amount of workers in a concurrent job
          Default is 1.
       num_cores : int, optional
          Number of cores which the GFN-XTB2 code
          will use.
       threshold: str, optional
          Level of geometry optimization required 
          as defined in GFN-XTB2 code. 

       Return
       ------
       None : 
          A relaxed structure as produced by the GFN-XTB2 
          code.
       """
       
       ending = self.ending
       gmts = self._seek_files(ending)

       
       with ThreadPoolExecutor(max_workers=int(max_work)) as executor:
           for i in tqdm(gmts, desc="Systems relaxed"):
               executor.submit(Htp.xtb_ff, xtb_path=path_xtb, xyz_geom=i,
                               init_dir=os.path.dirname(i),
                               num_cores=int(cores),
                               threshold=str(threshold))


    def htp_xtb_gfn(self, path_xtb, max_work, cores, threshold="crude"):
       """Function to perform a parallel concurrent calculation
           at the gfn-xtb level of theory.

       Parameters
       ----------
       path_xtb : str
          Path where the xtb executable is located. 
       max_work : int
          Number of maximum amount of workers in a concurrent job
          Default is 1.
       num_cores : int, optional
          Number of cores which the GFN-XTB2 code
          will use.
       threshold: str, optional
          Level of geometry optimization required 
          as defined in GFN-XTB2 code. 

       Return
       ------
       None : 
          A relaxed structure as produced by the GFN-XTB2 
          code.
       """
       
       ending = self.ending
       gmts = self._seek_files(ending)

       
       with ThreadPoolExecutor(max_workers=int(max_work)) as executor:
           for i in tqdm(gmts, desc="Systems relaxed"):
               executor.submit(Htp.xtb_gfn, xtb_path=path_xtb, xyz_geom=i,
                               init_dir=os.path.dirname(i),
                               num_cores=int(cores),
                               threshold=str(threshold))


    def htp_pyscf_semi(self, max_work, steps):
       """ Function to perform a parallel concurrent calculation
           at the gfn-ff level of theory.

       Parameters
       ----------
       max_work : int
          Number of maximum amount of workers in a concurrent job
          Default is 1.
       steps : int
          Number of steps which the PYSCF-semiempirical code
          will use to perform a geometry optimisation.

       Return
       ------
       None : 
          A relaxed structure as produced by the 
          PYSCF-semiempirical code.
       """
       
       ending = self.ending
       gmts = self._seek_files(ending)

       with ThreadPoolExecutor(max_workers=int(max_work)) as executor:
           for i in tqdm(gmts, desc="Systems relaxed"):
               executor.submit(Htp.pyscf_semi,steps, xyz_geom=i, save_path=os.path.dirname(i))


class Pyscf_print(object):
    """ Print a PYSCF object into cartesian 
        coordinates.
    """

    def xyz(self, mol, save_path):
        """Function to print a file with Cartesian 
           coordinates from a PYSCF geometry optimization
           calculation.

        Parameters
        ----------
        mol : pyscf.gto.mole.Mole
          A PYSCF gto mole Mole object.

        path_out : str
          Output path to print the relaxed coordinates.
       
        Result
        ------
        pysfc_final: str
          An xyz file containing the optimized
          coordinates of the provided PYSCF Mole
          object.
        """
        bohr2angs=0.529177
        coords = bohr2angs*mol.atom_coords()
        file_name="pyscf_final.xyz"
        complete_name=os.path.join(save_path, file_name)
        with open(complete_name, "w") as myfile:
             myfile.write('{}'.format(int(len(coords))))
             myfile.write('\n')
             myfile.write('\n')
             for i in range(len(mol.atom)):
                 myfile.write('{} {:.8f} {:.8f} {:.8f}\n'.format(str(mol.atom[i][0]),
                                                                 float(coords[i][0]),float(coords[i][1]),float(coords[i][2])))
             myfile.write('\n')
        print ("File pyscf_final.xyz has been created.")


