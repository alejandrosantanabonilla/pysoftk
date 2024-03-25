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
        print('Elapsed time for matrix calculation: {:.4f} seconds'.format(total_time))
        return result
    return timeit_wrapper


class micelle_whole(MDA_input):
   """A class used to make the largest micelle whole across 
      the pbc for better analysis.
   """

   def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """
       
      super().__init__(tpr_file, xtc_file)
  
   def obtain_largest_micelle_resids(self, spatial_clustering_results):
      """Function to obtain the largest micelle resids Bare in mind this may give two 
         different clusters if they have the same size - should check that you only select 
         one when using the results from here
         
         Note-need to check if it actually gives two different clusters if they have 
         the same size.
           
         Parameters
         ------------
         
         spatial_clustering_results : pandas.DataFrame
            Results from clustering analysis class. 
            Pandas dataframe with three columns: time (ps),  
            micelle_resids and micelle_size
          
             
         Returns
         --------
         
         longest_lists : list
            list with the number resids (of MDanalysis) of the largets micelle at the timesteps
            that want to be evaluated.
             
      """

      import pandas as pd
      import numpy as np

      
      df=pd.read_parquet(str(spatial_clustering_results))
      cluster_resids_tot = df['micelle_resids']
      longest_lists = [max(row, key=len) for row in cluster_resids_tot.values]

      return longest_lists



   def pbc_per_dimension(self, cluster_sel_positions, box, dimension, bins_min, bins_max, bins_step):
      """Function to find the broken parts of the micelle across the pbc and make it whole

         Parameters
         ------------
         
         spatial_clustering_results : np.array
            Numpy array with the positions of the resids 
            of interest. It needs to have the shape of n_molec x 3
         
         box : np.array
            Numpy array with the MDanalysis strat u.dimensions. Array with the x y and z dimensions
            follow by the angles of the box.

         dimension : class.int
            Dimension number that wants to be evaluated

         bins_min : class.int
            Smallest distance to evaluate if the pbc are broken (smallest x-y distance position)

         bins_max : class.int
            Biggest distance to evaluate if the pbc are broken (largest x-y distance position)

         bins_step : class.int
            Step between bins to evaluate broken pbc
            
         Returns
         --------
         
         cluster_sel_positions_return : np.array
            Array with the atom positions of the micelle made whole across the pbc
      
      """

      import MDAnalysis as mda 
      import MDAnalysis.analysis.distances
      import numpy as np
      import scipy.stats as stats

      cluster_sel_positions_return=cluster_sel_positions[:,dimension]

      y = stats.binned_statistic(cluster_sel_positions_return,
                                     cluster_sel_positions_return,
                                     bins=np.arange(bins_min, box[0]+bins_max, bins_step),
                                     statistic='count').statistic
        
       
      a = stats.binned_statistic(cluster_sel_positions_return,
                                   cluster_sel_positions_return,
                                   bins=np.arange(bins_min, box[0]+bins_max, bins_step),
                                   statistic='count').bin_edges
        

        
      move_above=False
      filled=np.where(y!=0)[0]

      move_above = [a[filled[i+1]] for i in range(len(filled)-1) if filled[i]!=filled[i+1]-1 ]

          
      if move_above!= False and len(move_above)>0:
         array = cluster_sel_positions_return

         for item in move_above:
            cluster_sel_positions_return[np.where(array > item)] = array[np.where(array > item)]-box[dimension] 

       
      return cluster_sel_positions_return


   def get_polymer_positions(self, u, frame, resname, cluster_resids):
      """Function to find the polymer positions

         Parameters
         ------------
         
         u : MDAnalysis.Universe
           An user-provided MDAnalysis universe.
         
         frame : class.int
            Time frame for calculation to be done.

         resname : class.list
            List of strings with the resname of the polymers.

         cluster_resids : class.str or class.int
            List with the polymer resids. 
            
             
         Returns
         --------
         
         box : np.array
            Array with the u.dimensions

         cluster_atom_positions : np.array
            Array with the atom positions of the micelle/cluster.

         cluster_sel : np.array
            Array with the atoms that belong to the micelle/cluster.
         
      """

      import MDAnalysis as mda
    
      u.trajectory[frame]
      #print(" ".join([str(item) for item in resname]))
      #print(" ".join([str(i) for i in cluster_resids]))
      cluster_sel = u.select_atoms('resname '+" ".join([str(item) for item in resname]) + ' and resid '+" ".join([str(i) for i in cluster_resids]))
      #print(cluster_sel) 
      cluster_atoms_positions=cluster_sel.positions.copy()
      box=u.dimensions  

      return box, cluster_atoms_positions, cluster_sel

   def count_dimensions(self, box, cluster_atoms_positions, dimensions):
      """Function to get the dimensions of a system, the box size and 
         cluster_atom_positions in the form to use in make_cluster_whole

         Parameters
         ------------
         
         box : np.array
            u.dimensions
        
         cluster_atom_positions : np.array
            Array with the atom positions of the polymers belonging to the cluster

         dimensions :  class.int
            Number of dimensions
             
         Returns
         --------

         box : np.array
            Array with the u.dimensions to feed in into make_cluster_whole

         cluster_atom_positions : np.array
            Array with the atom positions of the micelle/cluster in the correct shape 
            to feed it into make_cluster_whole

         dimensions : class.list
            List with the number of dimensions. E.g. 3D=[0, 1, 2]
         
      """
      dimensions = [i for i in range(dimensions)]
      cluster_atoms_positions = len(dimensions)*[cluster_atoms_positions]
      box = len(dimensions)*[box]

      return box, cluster_atoms_positions, dimensions

   def make_cluster_whole(self,
                       frame,
                       resname,
                       cluster_resids,
                       bins_min,
                       bins_max,
                       bins_step,
                       dimensions=3
                       ):
      
      """Function to make a polymer cluster whole across the pbc

         Parameters
         ------------
         
         frame : class.int
            Time frame that wants to be made whole

         resname : class.list
            List with the resname of the molecules of the cluster

         cluster_resids : class.list
            List with resids of the polymers that want to be made whole

         bins_min : class.int
            Smallest distance to evaluate if the pbc are broken (smallest x-y distance position)

         bins_max : class.int
            Biggest distance to evaluate if the pbc are broken (largest x-y distance position)

         bins_step : class.int
            Step between bins to evaluate broken pbc
            
         
         Returns
         --------
         
         frame : class.int
            Time frame where the micelle has been made whole

         atom_positions_whole : np.array
            Array with the atom positions of the micelle made whole across the pbc
      
      """
      import MDAnalysis as mda
      import numpy as np
      
      u=super().get_mda_universe()

      u.trajectory[frame] 
      box, cluster_atoms_positions, cluster_sel = self.get_polymer_positions(u, frame, resname, cluster_resids)
      res_box, res_clu_atoms, res_dimensions=self.count_dimensions(box, cluster_atoms_positions, dimensions)
      
      bins_min = len(res_dimensions)*[bins_min]
      bins_max = len(res_dimensions)*[bins_max]
      bins_step = len(res_dimensions)*[bins_step]
    
      atom_positions = list(map(self.pbc_per_dimension, res_clu_atoms, res_box, res_dimensions, bins_min, bins_max, bins_step))
      atom_positions_whole = np.stack(atom_positions, axis = 1) 
   
      return (frame, atom_positions_whole)


   @timeit
   def running_make_cluster_whole(self, resname, cluster_resids, start, stop, skip, bins_min=-50, 
                                  bins_max=50, bins_step=10):
      """Function to run the make_cluster_whole function over several frames 

         Parameters
         ------------
         
         resname : class.list
            List with the resname of the molecules of the cluster in each time step.

         cluster_resids : class.list
            List with resids of the polymers that want to be made whole in each time step.

         start : class.int
            Starting frame to perfomr the calculation.

         stop : class.int
            Last frame to perform the calculation.

         skip : class.int
            Length of step to perform the calculation every number 
            of frames.

         bins_min : class.int
            Smallest distance to evaluate if the pbc are broken (smallest x-y distance position).

         bins_max : class.int
            Biggest distance to evaluate if the pbc are broken (largest x-y distance position).

         bins_step : class.int
            Step between bins to evaluate broken pbc.
              
             
         Returns
         --------
         
         atom_positions_over_trajectory :  class.list 
            List of np.arrays where each entry of the list is the position of atoms that are 
            made whole across the pbc at a specific time step.

      """

      import MDAnalysis as mda
      import numpy as np
      
      from tqdm.auto import tqdm

      
      u=super().get_mda_universe()
      frames =[ts.frame for ts in u.trajectory[int(start):int(stop):int(skip)]]

      resname = [resname]*len(frames)
      cluster_resids_f = [cluster_resids]
                                      
      bins_min=len(frames)*[bins_min]
      bins_max=len(frames)*[bins_max]
      bins_step=len(frames)*[bins_step]

      # Single-threaded version of
      #atom_positions_over_trajectory = executor.map(self.make_cluster_whole, frames, resname, cluster_resids_f[0], bins_min, bins_max, bins_step)
                                      
      
      atom_positions_over_trajectory = tqdm(map(self.make_cluster_whole, frames, resname, cluster_resids_f[0], 
                                                                bins_min, bins_max, bins_step),total=len(frames))

      return list(atom_positions_over_trajectory)


   def get_atom_name_list(self, atom_sel):
      """Function to obtain the residues based on user-provided atom selection.

        Parameter
        -----------

        atom_sel: class.str
          User Provided name/tag of the molecules inside the box.

        Returns
        ---------

        final: class.list
          List containing all names of atoms.

      """
      import itertools
        
      b = [[str(item), str(item)] for item in atom_sel]
      c = sum(b, [])

      #creating permutation list
      final=list(set(list(itertools.combinations(c, 2)))) 

      return final


   def obtain_snapshot(self, output_name, atom_positions_at_specific_frame, cluster_resids, polymer_resname, 
                       frame, water_resname= ['SOL'], water_atom_name = ['OW'], solvant=False):
      """Function to obtain an snapshot of the micelle (pdb) file of a specific 
         frame and returns the new universe to analyse.
         
         Parameters
         -----------

         output_name : class.str
            Name of output file

         atom_positions_at_specific_frame : np.array
            Position of atoms at a specific frame

         cluster_resids : class.list
            List of polymers resids for the snapshot

         frame : class.int
            Time frame for snapshot

         water_resname : class.list
            List with name of solvants

         water_atom_name : class.list
            List of solvant atom names that want to be outputed in the sanpshot. E.g. with water 
            you may only be interested in outputing the oxygen water atom for further calculations

         solvant : class.boolean
            Set to True if interested in obtainig an output snapshot with the solvant too. 

         
         Returns
         --------

         u3: MDAnalysis.Universe
            An MDAnalysis universe containing the resulting frames

         output: MDAnalysis.Write
            A PDB file containing the selected molecules.
            
      """

      import MDAnalysis as mda
      import numpy as np
      import concurrent.futures
      from MDAnalysis import transformations as trans
  
      
      u2=super().get_mda_universe()

      u2.trajectory[int(frame)]
    
      positions_atoms = list(atom_positions_at_specific_frame[1])
      micelle = u2.select_atoms('resname '+" ".join([str(i) for i in polymer_resname])+' and resid '+" ".join([str(i) for i in cluster_resids]))

      print(len(micelle))
      print(len(positions_atoms))
     
      micelle.positions = positions_atoms


      #if any of the box is 

      if solvant == True:         
         system = u2.select_atoms('resid '+" ".join([str(i) for i in cluster_resids])+' or (resname '+" ".join([str(item) for item in water_resname]) +' and name '+" ".join([str(item) for item in water_atom_name])+')')
   
         workflow = (trans.unwrap(system),
                     trans.center_in_box(micelle, center='mass'),
                     trans.wrap(system, compound='fragments'))

         u2.trajectory.add_transformations(*workflow)
    
         with mda.Writer(str(output_name), system.n_atoms) as W: 
            W.write(system)

         u3 = mda.Universe(str(output_name))
         return u3
      

      else:
          with mda.Writer(str(output_name), micelle.n_atoms) as W:
               W.write(micelle)

          u3 = mda.Universe(str(output_name))
          return u3
