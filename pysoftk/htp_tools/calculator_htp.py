import glob
import subprocess as sp
from   subprocess import call
from concurrent.futures import ThreadPoolExecutor, as_completed
import os, fnmatch
from tqdm import tqdm
import time

class Htp(object):
    """Class enabling High-throughput calculations in different 
       created folders for provided molecules.

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

    This class requires tqdm, GFN-XTB and/or PYSCF semiempirical installation to be used.
 
    """

    def __init__(self, ending):
      """Parameters to initialize the HTP process

        Parameters
        ----------

        ending : str
           Name of the extension to be browsed in a provided folder. 

      """
      self.ending = ending


    @staticmethod   
    def xtb_gfn(xtb_path, xyz_geom, init_dir, num_cores=None, threshold=None):
       """Function invoking GFN-XTB Force-Field implementation for geometry optimization.

       Parameters
       -----------

       xtb_path : str
           Path where the GFN-XTB2 executable is located

       num_cores : int, optional
             Number of cores which the GFN-XTB2 code will use.

       threshold : str, optional
             Level of geometry optimization required as defined in 
             GFN-XTB2 code.

       xyz_geom : str 
             Provided geometry to be relaxed.

       init_dir : str
             Root path directory where folders containing the structures are found.
 
       Returns
       -------

       output : str
          File name output.log containing the optimization process carried out by GFN-XTB2.
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
       -----------

       xtb_path : str
           Path where the GFN-XTB2 executable is located

       num_cores : int, optional
           Number of cores used by the GFN-XTB code.

       threshold : str, optional
           Level of geometry optimization required as defined in GFN-XTB2 code.
 
       xyz_geom : str 
           Provided geometry to be relaxed.

       init_dir : str
           Root path directory where folders containing the structures are found.
 
       Returns
       --------

       output : str
           File name output.log containing the optimization process carried out by GFN-XTB2.
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

    def _find_files(self, directory, pattern):
        """Function to find files within a given directory based on a provided pattern.

        Parameters
        ----------

        directory : str
           Path of the main directory where the search will be carried out. 

        pattern : str
           Provided Alpha-pattern (extension of the molecular file) which will be used to perform the search.

        Results
        --------

        None :
            The results of the search as a string

        """
        for root, dirs, files in os.walk(directory):
            for basename in files:
                if fnmatch.fnmatch(basename, pattern):
                    filename = os.path.join(root, basename)
                    yield filename
       
    def _seek_files(self, file_ending):
        """Function to seek files within the current working directory based on a provided pattern.

        Parameters
        ----------

        file_ending : str
           Extension of the provided pattern to perform the search.

        Results
        -------

        gmts : list
           A list of the resulting files within the current working directory with the provided extension (ending). 

        """
        seek_file = "".join(("*",".", str(file_ending)))
        gmts = [filename for filename
                in self._find_files(os.getcwd(), str(seek_file))]
        return gmts


    def htp_xtb_ff(self, path_xtb, max_work, cores, threshold="crude"):
       """Function to perform a parallel concurrent calculation at the gfn-ff level of theory.

       Parameters
       ----------

       path_xtb : str
          Path where the xtb executable is located. 

       max_work : int
          Number of maximum amount of workers in a concurrent job default is 1.

       num_cores : int, optional
          Number of cores which the GFN-XTB2 code will use.

       threshold: str, optional
          Level of geometry optimization required as defined in GFN-XTB2 code. 

       Return
       ------

       None : 
          A relaxed structure as produced by the GFN-XTB2 code.

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
       """Function to perform a parallel concurrent calculation at the gfn-xtb level of theory.

       Parameters
       ----------

       path_xtb : str
          Path where the xtb executable is located. 

       max_work : int
          Number of maximum amount of workers in a concurrent job default is 1.

       num_cores : int, optional
          Number of cores which the GFN-XTB2 code will use.

       threshold: str, optional
          Level of geometry optimization required as defined in GFN-XTB2 code. 

       Return
       ------

       None : 
          A relaxed structure as produced by the GFN-XTB2 code.
       """
       
       ending = self.ending
       gmts = self._seek_files(ending)

       
       with ThreadPoolExecutor(max_workers=int(max_work)) as executor:
           for i in tqdm(gmts, desc="Systems relaxed"):
               executor.submit(Htp.xtb_gfn, xtb_path=path_xtb, xyz_geom=i,
                               init_dir=os.path.dirname(i),
                               num_cores=int(cores),
                               threshold=str(threshold))
