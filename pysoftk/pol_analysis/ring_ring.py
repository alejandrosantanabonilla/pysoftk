import os
import time
from functools import wraps

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.tools.utils_tools import *

# suppress some MDAnalysis warnings when writing PDB files
import warnings
warnings.filterwarnings('ignore')


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


class RSA(MDA_input):
    """A class used to compute a Ring Stacking Analysis (RSA).
    """
    
    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """

      super().__init__(tpr_file, xtc_file)

      
    def calc_angle(self, u_norm, u_norm2):
        """ Function to compute the dot product between two normal 
            vectors. 

        Parameters
        -----------

        u_norm: numpy.ndarray
           An user-provided normal vector. 

        u_norm2: numpy.ndarray
           An user-provided normal vector. 

        Return
        -------

        None:
           Result of a dot product between two vectors.

        """
        import numpy as np

        return abs(np.rad2deg(np.arccos(np.dot(u_norm,u_norm2)/np.linalg.norm(u_norm)/np.linalg.norm(u_norm2))))

    def GetRingSystems(self, min_1, max_1, min_2, max_2, includeSpiro=False):
        """Function to detect the atomic indexes  that  belong to a ring structure for a provided 
           molecular complex (single, dimer, timer and so on molecules).
 

           Parameters
           -----------

           min_1: int
              smallest atom index of molecule 1

           max_1: int
              largest atom index of molecule 1

           min_2: int
              smallest atom index of molecule 2

           max_2: int
              largest atom index of molecule 2

           includeSpiro: class.bool
              To check for spiro sites


           Returns
           --------
 
           ring_contact_indices: class.list
              A list of arrays with the mda atom indices of the 
              rings of the polymers.

        """
        import MDAnalysis as mda
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol=super().get_mda_universe()
        universe = mol.select_atoms('index '+str(int(min_1))+':'+str(int(max_1))+' or index '+str(int(min_2))+':'+str(int(max_2)))
        rdkit_mol = universe.convert_to("RDKIT")

        ri = rdkit_mol.GetRingInfo()
        
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon>1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
                    
            nSystems.append(ringAts)
            systems = nSystems 

        #getting the atom indices in mdanalysis
        ring_contact_indices = [universe.atoms[list(item)].indices for item in systems]    
      
        return ring_contact_indices 


    def separate_rings(self, ring_contact_indices, pol2_min):
        """Function to get the possible ring combinations between two polymers to 
           calculate the distances between them.
 
           Parameters
           -----------

           ring_contact_indices: class.list
              A list of arrays with the mda atom indices of the 
              rings of the polymers.

           pol2_min: int
              smallest atom index of the molecule with largest atom index 


           Returns
           --------
 
           ring_comb: class.list
              A list of arrays with the mda atom indices of the 
              rings of the polymers in contact.

        """
        import numpy as np
        import itertools

        pol1_indices = [item for item in ring_contact_indices if np.all(item < pol2_min)]
        pol2_indices = [item for item in ring_contact_indices if np.any(item >= pol2_min)]
      
        ring_comb = list(itertools.product(pol1_indices, pol2_indices))

        return ring_comb

    def get_polymer_resids(self, pol_indices):
        """Function to get the possible ring combinations between two polymers to 
           calculate the distances between them.
 
           Parameters
           -----------

           pol_indices: np.array
              Numpy array with the atom indices of the dimer.

           Returns
           --------
 
           ring_comb: class.list
             A list of arrays with the combination of mda molecule resids of the 
             two molecules involved

        """
        import itertools

        ring_comb = list(itertools.product(pol_indices[0], pol_indices[1]))
 
        return ring_comb 


    def svd(self, coordinates):
        """Function to obtain the Single Value Decomposition (SVD) 
           of coordinates.

        Parameters
        -----------
    
        coordinates: numpy.ndarray
           Array containing the coordinates upon which svd is going 
           to be calculated.
    
        Return
        -------
    
        None:
            A numpy.ndarray with the right unitary singular vector, orthogonal to 
            the given set of coordinates.

        """

        import numpy as np
        u_, s, vh = np.linalg.svd(coordinates - coordinates.sum(axis=0) / coordinates.shape[0])
   
        return vh[2, :] 

    
    def rings_stacking_dist(self, u, ring_comb, polymer_pos_1,
                            polymer_pos_2, pol_indices, cut_off_stack=5, par_scheme='OpenMP'):
        """Function to find if two rings are stacked. This is done by checking if 
           the distance of geometry between two rings is smaller than an 
           user-defined cut off distance.
 

            Parameters
            -----------

            ring_comb: numpy.ndarray
               Numpy array with all possible ring combinations

            polymer_pos_1: numpy.ndarray
               Numpy array with the correct positions of all atoms 
               of the polymer.

            polymer_pos_2: numpy.ndarray
                Numpy array with the correct positions of all atoms 
                of the polymer.

            pol_indices: numpy.ndarray
                Numpy array with the first atom index of each polymer. The 
                first one needs to be the one corresponding to polymer_pos_1
                and the second one to polymer_pos_2

            cut_off_stack : class.float
                Cut off distance for ring stackings


           Return
           --------
 
           vector_angle: class.float or class.boolean
                Angle between the two planes of atoms


        """
        import numpy as np
        import MDAnalysis as mda

        mol1_pos = polymer_pos_1[ring_comb[0]-pol_indices[0]]
        mol2_pos = polymer_pos_2[ring_comb[1]-pol_indices[1]]

        tot_dist = mda.lib.distances.distance_array(mol1_pos,
                                                    mol2_pos,
                                                    backend=str(par_scheme),
                                                    box=u.dimensions)

        
        if np.min(tot_dist) < cut_off_stack:      
            u_norm = self.svd(mol1_pos)
            u_norm2 = self.svd(mol2_pos)

            #finding the angle between the two vectors
            vector_angle = abs(np.rad2deg(np.arccos(np.dot(u_norm,u_norm2)/np.linalg.norm(u_norm)/np.linalg.norm(u_norm2))))

        else:
            vector_angle = None
            
        return vector_angle

    def run_rings_stacking_dist(self, u_frame, ring_comb, cut_off, ang_cutoff):
        """Function to find if two rings are stacked. This is done by checking 
           if the distance of geometry between two rings is smaller than a 
           user-defined cut off distance.
 

           Parameters
           -----------
           u_frame: int
               time step of the calculation

           ring_comb: numpy.ndarray
              Numpy array with all possible ring combinations

           cutoff : class.float
              Cut off distance for ring stackings

           ang_cut_off : class.float
              Cut off angle for ring stackings


        Returns
        --------
 
           ring_contact_atoms: class.list
              A list of arrays with the mda atom indices of the rings of the polymers

           pol_idx: class.list
               A list of the resids of the molecules with ring stacking

        """

        import MDAnalysis as mda
        import numpy as np
        from MDAnalysis.lib import distances 
        
        if len(ring_comb) < 1:            
            pass

        else:
            u=super().get_mda_universe()
            u.trajectory[u_frame]
            resid_1=np.unique(u.atoms[ring_comb[0][0]].resids)
            resid_2=np.unique(u.atoms[ring_comb[0][1]].resids)

            len_mol_1 = len(u.select_atoms('resid '+str(resid_1[0])))
            len_mol_2 = len(u.select_atoms('resid '+str(resid_2[0])))

            pol_resids = [list(resid_1), list(resid_2)]
            pol_resids_f = [item for sublist in pol_resids for item in sublist]
        
            #number to substract to the atom indices to only do the atom positions once
            pol1_first_index = u.select_atoms('resid '+str(pol_resids_f[0]))[0].index
            pol2_first_index = u.select_atoms('resid '+str(pol_resids_f[1]))[0].index
            pol_indices_f = [pol1_first_index, pol2_first_index]

            atom_pos_list = u.select_atoms('resid '+str(pol_resids_f[0])+'  '+str(pol_resids_f[1])).positions
            atom_positions_tot=atom_pos_list

            polymer_1_pos = atom_positions_tot[:len_mol_1]
            polymer_2_pos = atom_positions_tot[len_mol_1:]
    
            angles_tot_list= list(map(lambda i: self.rings_stacking_dist(u, ring_comb[i], polymer_1_pos,
                                                                polymer_2_pos, pol_indices_f, cut_off), range(len(ring_comb))))

            ring_index_cutoff = [i for i, x in enumerate(angles_tot_list)
                             if x is not None and (x < ang_cutoff or x > 180-ang_cutoff)]
        
            ring_contact_atoms = [ring_comb[i] for i in ring_index_cutoff]

            if len(ring_contact_atoms)>0:
               snapshot = u.select_atoms('resid '+str(pol_resids_f[0])+' '+str(pol_resids_f[1]))
               snapshot.positions = atom_positions_tot
            
               with mda.Writer('snapshot_'+str(pol_resids_f[0])+'_'+str(pol_resids_f[1])+'_'+str(u.trajectory.time)+'.pdb', snapshot.n_atoms) as W:              
                     W.write(snapshot)

               pol_idx = (pol_resids_f[0], pol_resids_f[1])

               return (ring_contact_atoms, pol_idx)


    def find_several_rings_stacked(self, rings_df_name):
        """Function to find several rings stacked using a provided
           pandas DataFrame from previous calculations.

        Parameters
        -----------

        rings_df_name: pandas.DataFrame
           Pandas dataframe containing the results from the ring_ring 
           analysis.


        Results
        --------

        connected_components: class.list
           List of connected polymers along the trajectory.

        """
        import networkx as nx
        import pandas as pd

        rings_df = pd.read_parquet(rings_df_name)
        
        connected_components_list=[]
        for i in range(len(rings_df)):
    
            flat_connection=[]   
            connections=rings_df.iloc[i, 1].tolist()
            flat_connection = [inner for outer in connections for inner in outer]
            graph = nx.Graph()

            # Add edges to the graph based on the connections
            for connection in flat_connection:
               elem1, elem2 = connection
               graph.add_edge(elem1, elem2)

            #Find the connected components in the graph
            connected_components = list(nx.connected_components(graph))
            connected_components_list.append(connected_components)

        return connected_components_list


    def nearest_neighbours(self, u, lab1, lab2):
        """Function to filter different molecules
           by closest interatomic distance.

        Parameters
        -----------

        u: MDAnalysis.object
           Universe containing the trajectory

        lab1: str
           User provided label of one molecule

        lab2: str
           User provided label for second molecule

        frame : int
            Frame for calculation

        Returns
        --------

        float: np.float
          Minimum distance between the atoms belonging to the polymer pair 
          and all the atompositions.

        """

        import MDAnalysis as mda
        import numpy as np
        from MDAnalysis.lib import distances 

        mol=u.select_atoms(str(lab1)) 
        mol1=u.select_atoms(str(lab2))

        tot_dist = distances.distance_array(mol.positions, mol1.positions, box=u.dimensions)

        return np.min(tot_dist) 

    
    def number_pol(self, u):
        """Number of polymers inside an user-provided box.

        Parameters
        -----------

        u: MDAnalysis.Universe
           An user-provided MDAnalysis universe.

    
        Return
        -------

        None:
          Number of polymers within a box.

        """

        import MDAnalysis as mda

        return len(u.atoms.residues) 

    
    def pol_cutoff(self, u, cutoff):
        """Function to count the number of polymers within an 
           user-provided cutoff.

        Parameters
        -----------

        u: MDAnalysis.universe
           An user-provided universe of MDAnalysis
 
        cutoff: float
           An user provided cutoff defining closest distances.

        frames : list
           User provided frames where calculation will be carried out

        Return
        -------

        index_cutoff: list
           List of intergers with the indexes of the selected
           molecules.

         b: list
            List of possible combinations of the selected molecules
     
        """
        import MDAnalysis as mda
        import numpy as np
        from tqdm import tqdm
        import itertools
        from itertools import combinations

        a=np.unique(u.atoms.resids)
        b=list(combinations(a, 2))  

        all_comb=[['{} {}'.format('resid',values[0]), '{} {}'.format('resid', values[1])]
                  for idx, values in enumerate(b)] 
      
        dist_tot=np.array([self.nearest_neighbours(u, str(i), str(j)) for i, j in tqdm(all_comb)])
      
        #returns list with the pair indices of the polymers that have at least two atoms that are in contact
        index_cutoff=np.where(dist_tot < float(cutoff))[0] 
    
        return index_cutoff,b
    
    @timeit
    def stacking_analysis(self, cut_off, ang_cut, start, stop, step, output_name="output.parquet"):
        """Function to perform a stacking analysis over several frames 
         Parameters
         -----------

         u: MDAnalysis.universe
           An user-provided universe of MDAnalysis
 
         cutoff: float
           An user provided cutoff defining closest distances.
         
         ang_cut: float
           An user provided cutoff defining the valid ring stacking angle range.

         start : int
           starting frame to perform the ring analysis.

         stop : int
           stopping frame of the ring analysis.

         step : int
           number of frames to skip during the ring analysis.

         output_name : nstr
            name of output parquet file with  the ring analysis.

         Return
         -------

          None: 
            pandas dataframe with atom indexes of rings where there is ring stacking at each time frame.
     
        """
   
        import pandas as pd
        import MDAnalysis as mda
        import numpy as np
        from tqdm import tqdm
        import itertools
        from itertools import combinations
        
        print('Ring Stacking analysis has started')       
        u=super().get_mda_universe()
        
        #Get the indices of b of the polymers that are in contact
        #Dimers with smallest distances according to radius cutoff
        index_cutoff,b=self.pol_cutoff(u, float(cut_off)) 
        dimer_mol=[b[i] for i in index_cutoff] #list of resnum of polymers that are in contact
        
        #Array with atom indices of the polymers that are in contact given in only one array
        dimer_ind=[u.select_atoms('resid {} {}'.format(values[0], values[1])).indices
                   for idx, values in enumerate(dimer_mol)]  

        dimers_indexes=[(min(u.select_atoms('resid {}'.format(values[0])).indices), 
        max(u.select_atoms('resid {}'.format(values[0])).indices), 
        min(u.select_atoms('resid {}'.format(values[1])).indices), max(u.select_atoms('resid {}'.format(values[1])).indices))
                   for idx, values in enumerate(dimer_mol)] 
       
        us = [u]*len(dimers_indexes)  

        dim2 = list(tqdm(map(self.GetRingSystems, [x[0] for x in dimers_indexes],
                             [x[1] for x in dimers_indexes], [x[2] for x in dimers_indexes],
                             [x[3] for x in dimers_indexes]), total=len(dimers_indexes),
                             desc="Detecting atoms in Rings"))
    
        u4 = [u]*len(dim2)

        dim3 = tqdm(map(self.separate_rings, dim2, [x[2] for x in dimers_indexes]),
                    total=len(dim2), desc="Separating Rings")

        c = list(dim3)
        u_f=[ts.frame for ts in u.trajectory[int(start):int(stop):int(step)]]

        print ("\n")
        print ("Preparing DataFrame to store the results")
          
        data = [[[], []] for _ in range(len(u_f))]
        timestep_count=0

        for u_t in tqdm(u_f):           
            u4 = [u_t]*len(c)
            cuts = [float(cut_off)]*len(c)
            ang_cuts = [float(ang_cut)]*len(c)

            dim4 = list(tqdm(map(self.run_rings_stacking_dist, u4, c, cuts, ang_cuts),
                         total=len(dim3), desc="Computing stacking distances"))
    
            results = [result for result in dim4 if result is not None]

            if len(results)<1:
                  print('No ring stacking in this frame')
            else:
               results1, results2 = zip(*results)
               print(results1)
                
               flattened_list_1 = []
               for item in results1:
                     for array in item:
                           flattened_list_1.extend(array)
             
              
               data[timestep_count][0].append(flattened_list_1)
               data[timestep_count][1].append(results2)
            timestep_count+=1
            
        df_results = pd.DataFrame(data, columns=['atom_index', 'pol_resid'])

        # Save the DataFrame to a file
        df_results.to_parquet(str(output_name), index=False)
    
        print ("Information succesfully stored in {}".format(str(output_name)))
        print ("Stacking analysis has succesfully finished!")
    
