import os
import time
from functools import wraps
from pysoftk.pol_analysis.tools.utils_mda import MDA_input
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

    
class SCP(MDA_input):
    """
    A class used to compute a spatial cluster for polymers (SCP)
    """
    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
      
      super().__init__(tpr_file, xtc_file)
      

    def contact_matrix(self, u, name_list, cluster_cutoff):
        """ Function to compute the contact matrix.

        Parameters:
        ------------

        u: MDAnalaysis.Universe
           An user-provided MDAnalysis universe.

        name_list: class.list
           A list containing the names used for the contact 
           matrix analysis.

        cluster_cutoff: class.float
           A user-defined number for a radius cutoff.


        Returns:
        ----------

        pairs: class.list
           A list of tuples containing the indexes of the selected 
           molecules.
           
        """
        import MDAnalysis as mda 
        import numpy as np

        atoms_group1=u.select_atoms('name '+str(name_list[0])).atoms
        atoms_group2=u.select_atoms('name '+str(name_list[1])).atoms
        group1=u.select_atoms('name '+str(name_list[0])).atoms.positions
        group2=u.select_atoms('name '+str(name_list[1])).atoms.positions


        distances_ar = distance_matrix_calc(u, group1, group2)
        dist_array = np.full((len(u.select_atoms('name '+str(name_list[0]))), len(u.select_atoms('name '+str(name_list[1])))), False)
        dist_array[np.where(distances_ar <= cluster_cutoff)] = True
    
        source, target = np.where(dist_array)
        
        #pairs = [[s+1, t+1] for s, t in zip(source, target) if s != t]
        
        pairs = [(atoms_group1[source[i]].resid, atoms_group2[target[i]].resid) for i in range(len(source)) if source[i] != target[i]]

    
        return pairs


    def spatial_clustering(self, u, atom_sel, cluster_cutoff, frame):
        """Function that gives you the resids and sized of each 
           polymer aggregate per frame in a trajectory.

           Parameters
           -----------

           frame: MDanalysis.Object
              An MDAnalysis object that contains one frame of a trajectory.

           atom_sel: class.str
              User-provided string containing the ResID/tag of the molecule in the trajectory.  

           cluster_cutoff: class.int
              User defined cutoff to select as neighbour.

           max_workers: class.int
              Number of threads to be used in the calculation.

           Returns
           --------
   
           tuple containing clusters and cluster_sizes.

           clusters: class.list
               A tuple containing the ResIDs. 

           cluster_sizes: class.list
              Length of the computed clusters.
 
        """
        
        import MDAnalysis as mda

        u=u[0]

        u.trajectory[frame]

        name_list = get_atom_name_list(atom_sel)
        size_name_list = [cluster_cutoff]*len(name_list)

        func=lambda y,z: self.contact_matrix(u, y, z)
        pairs = map(func, name_list, size_name_list)

        pairs_flatten = sum(list(pairs), [])
        list_of_tuples = [tuple(lst) for lst in pairs_flatten]
   
        G_net=create_network(list_of_tuples)
        clusters=[h for h in G_net]
        cluster_sizes=[len(h) for h in clusters]

        # Providing an occurrence list based on cluster_sizes length
        no_ones=len(u.select_atoms('name '+str(name_list[0][0])))-sum(cluster_sizes)
        [cluster_sizes.append(1) for i in range(no_ones)] 

        return u.trajectory[frame].time, clusters, cluster_sizes

    @timeit
    def spatial_clustering_run(self, start, stop, skip, atom_names, cluster_cutoff, name_file):
       """Function to compute the spatial clustering function over an MDAnalysis
          trajectory. 

        Parameters
        ============

        start : class.int
            Frame to start calculation.

        stop : class.int
            Frame to stop calculation.
        
        step : class.int
            Number of frames to skip.

        atom_names: class.str
            Atomic names used to procure the search.

        cluster_cutoff: class.float
            User-provided value to define the cluster length.

        name_file: class.str
            Name of the User-supplied file.

        Returns
        --------
 
        None : 
            A pandas.DataFrame.parquet file with the results stored.

       """
       import MDAnalysis as mda
       import pandas as pd
       from tqdm.auto import tqdm
       
       u=super().get_mda_universe()
       frames = [ts.frame for ts in u.trajectory[int(start):int(stop):int(skip)]]

       u_f=get_universe_per_frame(u, start, stop, skip)

  
       func=lambda z: self.spatial_clustering(u_f, atom_names, cluster_cutoff, z)
       results = list(tqdm(map(func, frames), total=len(frames)))

       df = pd.DataFrame(results, columns=['time', 'micelle_resids', 'micelle_size'])
       final_name=[str(name_file),'.parquet']
       df.to_parquet(''.join(final_name))
       
       print ("The file {} has been successfully created.".format(''.join(final_name)))
