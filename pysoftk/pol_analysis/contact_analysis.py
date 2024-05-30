import os
import time
from functools import wraps

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *

from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
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

class contacts(MDA_input):
    """ A class used to compute the contacts between the polymers 
        of a micelle. 
    """

    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
      
      super().__init__(tpr_file, xtc_file)


    def contacts_calc(self, u, micelle_selection, micelle_positions, group1, group2, contact_distance):
        """Function to get atom pairs
 
           Parameters
           -----------

           u: MDAnalysis.Universe
               An user-provided MDanalysis trajectory.

           micelle_selection: MDAnalysis.Object 
                MDAnalysis atom group of atoms belonging to the micelle
            
           micelle_positions : np.array
                positions of all atoms from micelle selection

           group1 : class.str
                atom  names of the atoms for group1 of the contacts calculation

           group2 : class.str
                atom names of the atoms for group2 of the contacts calculation

           Returns
           --------
 
           contacts_mat: np.array
                Array with the contacts between polymers

        """

        import numpy as np
        import MDAnalysis as mda


        micelle = micelle_selection
        micelle.positions = micelle_positions[:][1]

        molec_1=micelle.select_atoms('name '+str(' '.join(group1)))
        molec_2=micelle.select_atoms('name '+str(' '.join(group2)))
        
        molec1_positions = molec_1.positions
        molec2_positions = molec_2.positions

        resids_1 = np.unique(molec_1.resids)
        resids_2 = np.unique(molec_2.resids)

        
        distance_matrix = distance_matrix_calc(u, molec1_positions, molec2_positions)

        indices_contact = np.where(distance_matrix <= contact_distance)
        indices_non_contact = np.where(distance_matrix > contact_distance)


        distance_matrix[indices_contact] = 1
        distance_matrix[indices_non_contact] = 0


        len1 = int(len(molec_1)/len(resids_1))
        len2 = int(len(molec_2)/len(resids_2))


        pairs = get_index_pairs(len(distance_matrix),
                                len1, len(distance_matrix[0]),
                                len2)

        distance_matrix_f = [distance_matrix]*len(pairs)
        len1_f = [len1]*len(pairs)
        len2_f = [len2]*len(pairs)
    
        mat = map(find_values, distance_matrix_f, pairs, len1_f, len2_f)

        contacts_mat = np.sum(list(mat), axis=0)


        return contacts_mat


    def run_contacts_calc(self,  micelle_selection, micelle_positions, group1, group2, contacts_distance):
        """Function to run contacts_calc function over timme
 
            Parameters
            -----------

            micelle_selection: class.list
                List with resids of the polymers belonging to the micelle

            micelle_positions : np.array
                Array with the positions of all atoms of the micelle

            group1 : class.str
                Atom  names of the atoms for group1 of the contacts calculation

            group2 : class.str
                Atom names of the atoms for group2 of the contacts calculation

            contacts_distance: class.int
                Cut off distance for contacts

            Returns
            --------
 
            contacts_over time: np.array
                Array with the contacts between polymers in each timestep

        """
        import MDAnalysis as mda
        import concurrent.futures
        from tqdm.auto import tqdm

        u=super().get_mda_universe()
        
        u_sel = [u]*len(micelle_positions)
        micelle_selection_f = map(selecting_atoms, u_sel, micelle_selection)

        group1_f=len(micelle_positions)*[group1]
        group2_f=len(micelle_positions)*[group2]

        contacts_distance_f = len(micelle_positions)*[contacts_distance]
        contacts_over_time = list(tqdm(map(self.contacts_calc, u_sel, micelle_selection_f, micelle_positions, group1_f, group2_f, contacts_distance_f), total=len(u_sel)))

        return contacts_over_time
  
