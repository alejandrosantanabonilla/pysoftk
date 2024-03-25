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

class intrinsic_density_water(MDA_input):
    """ A class used to compute the contacts between the polymers 
        of a micelle. 
    """

    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
      
      super().__init__(tpr_file, xtc_file)


    def select_water_atoms(self, u, name):

        """Function to select the solvent atoms

            Parameters
            ----------

            u: MDAnalysis.Universe
               An user-provided MDanalysis trajectory

            name : list
                List with the atom names of the solvent
        """


        water_selections = u.select_atoms('name '+str(' '.join(name)))

        return water_selections


    def intrinsic_density_not_normalized(self, u, frame, whole_micelle_largest_cluster_resids, micelle_positions, water_selection, core_group, n_rand_points, n_bins, n_min, n_max, n_step ):
        """Function to run the intrinsic density in a specific frame
 
           Parameters
           -----------

           u: MDAnalysis.Universe
               An user-provided MDanalysis trajectory.

           frame : int
                trajectory frame

           whole_micelle_largest_cluster_resids: list
                list with resids of polymers forming the desired micelle

            
           micelle_positions : np.array
                positions of all atoms from micelle selection

           water_selection: list
                list of solvent atoms to calculate the density

            core_group : class.str
                atom  names of the atoms of the core for density calculation

            n_rand_points : int
                number of water molecules

            n_bins : int 
                nnumber of bins for the density calculation

            n_min : int
                min bin value for density calculation

            n_max : int
                max bins value for density calculation

            n_step : int
                size of bin


           Returns
           --------
 
           intrinsic_r_total : np.array
                not normalized intrinsic density

           intrinsic_r_total_norm : np.array
                normalized intrinsic density

        """

        import numpy as np
        import MDAnalysis as mda
        import pysoftk.pysoftwhere.icsi as pysw_icsi
        import scipy.stats as stats
        from MDAnalysis import transformations as trans

        u2=super().get_mda_universe()
        u2.trajectory[frame]

        micelle = selecting_atoms(u2, whole_micelle_largest_cluster_resids)
        micelle.positions = micelle_positions[:][1]

        water=self.select_water_atoms(u2, water_selection)      
        system=water+micelle

        workflow=(trans.unwrap(system),
                   trans.center_in_box(micelle, center='mass'),
                   trans.wrap(system, compound='fragments'))

        u2.trajectory.add_transformations(*workflow)

        core=u2.select_atoms('name '+str(' '.join(core_group)))

        shell=water
        core_positions = core.positions
        shell_positions = shell.positions
        cluster_atoms_positions=micelle.positions
        
        intrinsic_r, spherical_r, icsi_vals = pysw_icsi.icsi(u,
                                    cluster_resids=whole_micelle_largest_cluster_resids,
                                    cluster_atoms_positions=cluster_atoms_positions,
                                    core_sel_atoms_positions=core_positions,
                                    shell_sel_atoms_positions=shell_positions,
                                    frame=frame,
                                    no_bins=n_bins,
                                    no_random_points=n_rand_points,
                                    normalisation_run=False)()
    
        intrinsic_r_total=stats.binned_statistic(intrinsic_r,intrinsic_r,bins=np.arange(n_min,n_max,n_step),statistic='count').statistic
        

        intrinsic_r_norm, spherical_r_norm, icsi_vals_norm = pysw_icsi.icsi(u,

                                    cluster_resids=whole_micelle_largest_cluster_resids,
                                    cluster_atoms_positions=cluster_atoms_positions,
                                    core_sel_atoms_positions=core_positions,
                                    shell_sel_atoms_positions=shell_positions,
                                    frame=frame,
                                    no_bins=n_bins,
                                    no_random_points=n_rand_points,
                                    normalisation_run=True)()
    
        intrinsic_r_total_norm=stats.binned_statistic(intrinsic_r_norm, intrinsic_r_norm, bins=np.arange(n_min,n_max,n_step),statistic='count').statistic


        return (intrinsic_r_total, intrinsic_r_total_norm)


    def box_volume(self, u, frame):

        """ Function to calculate the volume of the system box at a specific frame
        
            Parameters
            -----------

            u: MDAnalysis.Universe
               An user-provided MDanalysis trajectory.

            frame : int
                frame selected for volume calculation


            Returns
            --------
            vol : float
                volume of system box

            """

        u.trajectory[frame]

        vol=u.dimensions[0]*u.dimensions[1]*u.dimensions[2]

        return vol


    def run_intrinsic_density(self,  micelle_selection, micelle_positions, water_selection, core, heavy_atom_water_name, start, stop, step, n_bins=31, n_min=-40.5, n_max=150, n_step=0.1):
        """Function to run the intrinsic density over time
 
           Parameters
           -----------

           micelle_selection: list
                list of atoms belonging to the micelle
            
           micelle_positions : np.array
                positions of all atoms from micelle selection

            water_selection : list
                solvent atom selection

            core : class.str
                atom  names of the atoms of the core for density calculation

           shell : class.str
                atom names of the atoms of the shell for density calculation

            heavy_atom_water_name : class.str
                name of heavy water atoms

            start : int
                starting frame to perform calculation

            stop : int
                last frame to perform calculation

            step : int
                frames skipped in calculation

            n_bins : int 
                nnumber of bins for the density calculation

            n_min : int
                min bin value for density calculation

            n_max : int
                max bins value for density calculation

            n_step : int
                size of bin


           Returns
           --------
 
           final_density : np.array
                intrinsic density

            binned_space : np.array
            array with the values of the binned radial distance, this allows easier plotting of the density

        """
        import MDAnalysis as mda
        import concurrent.futures
        from tqdm.auto import tqdm
        import  numpy as np

        u=super().get_mda_universe()

        
        u_sel = [u]*len(micelle_positions)
        
        whole_micelle_largest_cluster_resids=micelle_selection
            
        n_bins_f=len(micelle_positions)*[n_bins]
        n_min_f=len(micelle_positions)*[n_min]
        n_max_f=len(micelle_positions)*[n_max]
        n_step_f=len(micelle_positions)*[n_step]
        
        core_f=len(micelle_positions)*[core]
      
        water=u.select_atoms('name '+str(' '.join(heavy_atom_water_name)))
        
        n_rand_points = 3*len(water)
        n_rand_points_f = len(micelle_positions)*[n_rand_points]
        frames=get_frames_hydr(u, start, stop, step)
        box_vol=list(map(self.box_volume, u_sel, frames))
        water_selection_f=[water_selection]*len(micelle_positions)
   
        intrinsic_density_not_normalized_over_time = list(tqdm(map(self.intrinsic_density_not_normalized, u_sel, frames, whole_micelle_largest_cluster_resids, micelle_positions, water_selection_f, core_f, n_rand_points_f, n_bins_f, n_min_f, n_max_f, n_step_f), total=len(u_sel)))

        id_not_normalized = [item[0] for item in intrinsic_density_not_normalized_over_time]
        id_normalized = [item[1] for item in intrinsic_density_not_normalized_over_time]

        intrinsic_density_profile=np.sum(np.array(id_not_normalized), axis=0)/len(micelle_positions)
        S_bar=np.sum(np.array(id_normalized), axis=0)*np.mean(box_vol)/(len(micelle_positions)*n_rand_points)
        final_density=intrinsic_density_profile/S_bar
        binned_space_plot=np.arange(n_min, n_max, n_step)[:-1]

        return final_density, binned_space_plot
