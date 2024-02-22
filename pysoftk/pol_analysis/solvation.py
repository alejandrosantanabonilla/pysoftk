import os
import time
from functools import wraps

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP

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


class solvation(MDA_input):
    """A class used to compute the contacts between the polymers of a micelle 
    """


    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
      
      super().__init__(tpr_file, xtc_file)


    def hydration_calc(self, frame, largest_cluster_resids,
                       micelle_pos, water_name, polymer_oxygen_names, cut_off):
        """Function to calculate the solvation of selected polymer atoms

        Parameters
        -----------

        frame : class.int
            Time frame of step for calculation.

        largest_cluster_resids : class.list
            List with resids of polymers forming the micelle at a 
            specific time step.

        micelle_pos : np.array
            3D np.array with the positions of the atoms that conform 
            the micelle whole (not broken across the pbc boundaries).

        water_name: class.str
            Water atom names for hydration calculation.

        polymer_oxygen_names : class.list
            List of atom names of the polymer for hydration calculation.

        cut_off : class.float
            Number cut off of the distance between the polymer beads 
            and the water for hydration calculation.

        Returns
        --------
 
        coord_number: list
            Coordination number of the polymer atoms.

        """

        import MDAnalysis as mda
        from MDAnalysis import transformations as trans

        u2 = super().get_mda_universe()

        #setting the universe in the right time step
        u2.trajectory[frame]
        micelle = selecting_atoms(u2, largest_cluster_resids)

        micelle.positions = micelle_pos[:][1]

        water=u2.select_atoms('name '+str(' '.join(water_name)))
        system = micelle + water

        #MDAnalysis workflow to compute a series of functions in one sequence.
        workflow = (trans.unwrap(system),
                   trans.center_in_box(micelle, center='mass'),
                   trans.wrap(system, compound='fragments'))
    
        u2.trajectory.add_transformations(*workflow)
        polymer_oxygens=micelle.select_atoms('name '+str(' '.join(polymer_oxygen_names)))
        polymer_oxygens_pos = polymer_oxygens.positions

        coord_number = coord_matrix_calc(u2, polymer_oxygen_names, polymer_oxygens_pos,
                                         water.positions, cut_off)
    
        return coord_number

    def solvation_calc_run(self, start_frame, stop_frame, step_frame, largest_cluster_resids, micelle_pos, 
                           water_name, polymer_oxygen_names, cut_off):
        """Function to calculate the solvation calculation over the selected frames 

        Parameters
        -----------

        start_frame : class.int
            frame to start calculation

        stop_frame : class.int
            frame to stop calculation
        
        step_frame : class.int
            number of frames to skip

        largest_cluster_resids : class.list
            list with resids of polymers forming the micelle at a specific time step

        micelle_pos : np.array
            3D np.array with the positions of the atoms that conform the micelle whole 
            (not broken across the pbc boundaries)

        water_name: class.str
            warer atom names for hydration calculation

        polymer_oxygen_names : class.list
            list of atom names of the polymer for hydration calculation

        cut_off : class.float
            number cut off of the distance between the polymer beads and the water for 
            hydration calculation.

        Returns
        --------
 
        coord_number: class.list
            coordination number of the polymer atoms at a specific frame

        """
        
        import numpy as np
        import MDAnalysis as mda 
        from tqdm.auto import tqdm

        u=super().get_mda_universe()
        frames = get_frames_hydr(u, start_frame, stop_frame, step_frame)
    
        water_name_f = [water_name]*len(frames)
        polymers_oxygen_names_f = [polymer_oxygen_names]*len(frames)
        cut_off_f = [cut_off]*len(frames)
        coord_number=list(tqdm(map(self.hydration_calc, frames, largest_cluster_resids, micelle_pos, water_name_f,
                                            polymers_oxygen_names_f, cut_off_f), total=len(frames)))

        #It is returning the coord number in each time step for more flexibility
        return coord_number
