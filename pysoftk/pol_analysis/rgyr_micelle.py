import os
import time
from functools import wraps

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.tools.utils_tools import *

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        import time
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper


class rgyr(MDA_input):
    """A class used to compute the radius of gyration of micelles 
    """

    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
      
      super().__init__(tpr_file, xtc_file)
      
    
    def rgyr_calc(self, atom_sel_micelle, micelle_positions, subgroup=[]):
        """Funcion to calculate the radius of gyration of the chosen atoms
 
        Parameters
        -----------

        atom_sel_micelle : MDAnalysis.atom_selection
             MDAnalysis atom_selection function.

        subgroup: MDAnalysis.atom_group
            MDAnalysis atom group to calculate the rgyr
            
        micelle_positions : np.array
            np.array with the positions of the micelle atoms

        Returns
        --------
 
        radius_of_gyration: class.float
            Radius of gyration of the selected atoms

        """

        import MDAnalysis as mda

        micelle = atom_sel_micelle
        micelle.positions = micelle_positions[:][1]

        if len(subgroup)>1:
            micelle_subgroup=micelle.select_atoms('name '+str(' '.join(subgroup)))

        else: 
            micelle_subgroup=micelle
   
        radius_of_gyration = (micelle_subgroup.radius_of_gyration())  

        return radius_of_gyration
    

    def running_rgyr(self, atom_names, micelle_positions,
                     start, stop, skip, subgroup=[]):
        """Function to calculate the radius of gyration of the 
           chosen atoms over time.

        Parameters
        -----------

        u: MDAnalysis.Universe
            An user-provided MDAnalysis universe.

        atom_sel_type : class.str
            Type of atom selection.

        atom_names : class.list
            List with the names of atoms.
            
        micelle_positions : np.array
            An array with the positions of the micelle atoms.

        subgroup: class.str
            A list of atom names to calculate the rgyr.
            
        Start : class.int
            Starting frame of the trajectory.

        Stop : class.int
            Ending rame of the trajectory.

        Skip : class.int
            Skipping every that many frames of the trajectory.

       Returns
       --------
 
       rgyr: np.array
           Radius of gyration of the selected atoms over time

       """

        import MDAnalysis as mda
        from tqdm.auto import tqdm
        
        u=super().get_mda_universe()
        subgroup_f = [subgroup]*len(u.trajectory[start:stop:skip])
     
        us = [u]*len(u.trajectory[start:stop:skip])         
        atom_sel_micelle = map(selecting_atoms, us, atom_names)

  
        rgyr_f = list(tqdm(map(self.rgyr_calc, atom_sel_micelle, micelle_positions, subgroup_f), total=len(us)))

        return rgyr_f
      
