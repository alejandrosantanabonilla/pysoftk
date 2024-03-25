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

class spherical_density_water(MDA_input):
    """ A class used to compute the contacts between the polymers 
        of a micelle. 
    """

    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
      
      super().__init__(tpr_file, xtc_file)


    def dist_com(self, u, frame, atom_sel_micelle, micelle_positions, water_selection):
    
    
        """Function to calculate the distance between atoms and the center of mass of the micelle
 

       Parameters
       -----------

       u: mda.Universe
            mda.Universe

       atom_sel_micelle : mda.object
            atom group selection of the micelle

        micelle_positions : np.array
            array with the positions of the micelle whole
            

        water_selection : mda.object
            atom group selection of the water


        

       Returns
       --------
 
       distance_com: np.array
            array with the distances of the aotms to the COM of the micelle

        """

        import MDAnalysis as mda
        import numpy as np
        import scipy.stats as stats
        from scipy.spatial.distance import cdist
        import MDAnalysis.analysis.distances

        



        u.trajectory[frame]
        micelle = atom_sel_micelle
        micelle.positions = micelle_positions[:][1] 

        micelle_com = np.array([micelle.center_of_mass()])

        atom_group_water = water_selection
    

        #distance_COM = cdist(atom_group.positions, micelle_com, 'euclidean')

        distance_COM=MDAnalysis.analysis.distances.distance_array(atom_group_water.positions, micelle_com, box=u.dimensions)

        
        return distance_COM


    def sph_histogram(self, distance, maxbin, minbin, step):

        """Function to create a histogram of the distances
 

       Parameters
       -----------

        distance: np.array
            distances to be binned

        maxbin : int
            largest bin value

        minbin : int
            smalles bin value

        step : int
            bin size



        

       Returns
       --------
 
       all_histograms: np.array
            array with the distances values binned

        """

        import numpy as np
        import scipy.stats as stats 
    
        x=len(np.arange(minbin,maxbin,step)[:-1])
    

        all_histograms=np.zeros(x)



        binned_space=np.arange(minbin,maxbin,step)
    
        binned_space=np.array(binned_space)


        l = [item for item in distance]

        flat_list = [item for sublist in l for item in sublist]

    
        all_histograms+=stats.binned_statistic(flat_list, flat_list, 'count', bins=binned_space).statistic
   
    
        all_histograms=np.array(all_histograms)


        return all_histograms



    def density_calc(self, binned_space, all_histograms):

        """Function to calculate the density
 

       Parameters
       -----------

        binned_space: np.array
            bins

       all_histograms : np.array
            distances histograms
      


        

       Returns
       --------
 
       density : np.array
            array with density

        """

        import numpy as np

        volume = 4/3 * (3.14*((binned_space[1]**3)-(binned_space[0]**3)))
        density = all_histograms/volume 
    
        return density


    def run_density_calc(self, micelle_sel_type, whole_micelle_selection, micelle_positions, water_sel_type, water_atom_selection, start, stop, step_frame, maxbin=200, minbin=0, step=0.1):


        """Function to get the mean  density calculation over time as a funciton wrt the micelle COM
 

       Parameters
       -----------

        micelle_sel_type : str
            string type for the selection of the micelle atoms

        whole_micelle_selection : np.array
            elements to select the micelle atoms

        micelle_positions: np.array
            array with the positions of the micelle atoms whole

        water_sel_type : str
            string type for the selection of the solvent atoms

        water_atom_selection : list
            list with the names of the solvent  chosen for the density calculation

        maxbin : int
            largest bin value

        minbin : int
            smalles bin value

        step : int
            bin size



        

       Returns
       --------
 
       density: np.array
            array with the density values as a function of the distance to the COM of the micelle
            
        binned_space : np.array
            array with the values of the binned radial distance, this allows easier plotting of the density

        """

        import numpy as np
        import matplotlib.pyplot as plt


        u=super().get_mda_universe()
    
        us = [u]*len(whole_micelle_selection)
        atom_sel_type_f = [micelle_sel_type]*len(whole_micelle_selection)
        water_sel_type_f = [water_sel_type]*len(whole_micelle_selection)
        water_atom_selection_f=[water_atom_selection]*len(micelle_positions)

        

        whole_micelle_selection_f = map(selecting_atoms_density, us, atom_sel_type_f, whole_micelle_selection)
        water_selection_f = map(selecting_atoms_density, us, water_sel_type_f, water_atom_selection_f)

        
        frames=get_frames_hydr(u, start, stop, step_frame)

       

        distance_COM_total = map(self.dist_com, us, frames, whole_micelle_selection_f, micelle_positions, water_selection_f)
        
        maxbin_f=len(micelle_positions)*[maxbin]
        minbin_f=len(micelle_positions)*[minbin]
        step_f=len(micelle_positions)*[step]
        binned_space=np.arange(minbin,maxbin,step)

       

        all_histograms = list(map(self.sph_histogram, distance_COM_total, maxbin_f, minbin_f, step_f))

    
        all_histograms_f = np.sum(all_histograms, axis=0)


        normalisation = len(micelle_positions)



        binned_space_f = create_binned_space(binned_space)

        density = map(self.density_calc, binned_space_f, all_histograms_f/normalisation)


        density_plot = list(density)

       


        return density_plot, binned_space[:-1]



   
