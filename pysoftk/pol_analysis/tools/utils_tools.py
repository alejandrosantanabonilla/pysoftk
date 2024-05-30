def selecting_atoms(u, atom_resids):
        """Funcion to calculate the radius of gyration of the chosen atoms
 
        Parameters
        -----------

         u: MDAnalysis.Universe
            An User-provided Universe.

        atom_sel_type : class.str
            Atom selection of the atoms that form the micelle for 
            rgyr calculation.
            
        atom_resids : class.list
            list with the resids of the polymers that form the micelle 
            for rgyr calculation.

        Returns
        --------
 
        selected_atoms: MDAnalysis.selection
            Atom selection

        """

        import MDAnalysis as mda
         
        selected_atoms = u.select_atoms('resid '+" ".join([str(i) for i in atom_resids]))

        return selected_atoms



def selecting_atoms_density(u, atom_sel_type, atom_names):

        """Function to select specific atoms
 

       Parameters
       -----------

       u: mda.Universe
            mda.Universe

       atom_sel_type : str
            type of atom selection, resid, names etc...
            

        atom_names : list
            list with the names of the atom selection


        

       Returns
       --------
 
       selected_atoms: mda.selection
            atom selection

        """

        import MDAnalysis as mda
     
        selected_atoms = u.select_atoms(str(atom_sel_type)+ ' '+" ".join([str(i) for i in atom_names]))


        return selected_atoms



def get_frames(tpr_file, xtc_file, start, stop, step):
        """Function to obtain the frames of a trajectory

        Parameters
        -----------

        tpr_file : MDAnalysis.tpr_file 
          A tpr trajectory to analyse.

        xtc_file : MDAnalysis.xtc_file 
          A XTC trajectory to analyse

        start : class.int
            frame to start calculation

        stop : class.int
            frame to stop calculation
        
        step : class.int
            number of frames to skip


        Returns
        --------
 
        frames : np.array
            np.array with the selected frames
           
        """

        import numpy as np
        import MDAnalysis as mda

        u=mda.Universe(tpr_file, xtc_file)

        return np.array(list(map(lambda ts: ts.frame, u.trajectory[start:stop:step])))


def get_frames_hydr(u, start, stop, step):
        """Function to obtain the frames of a trajectory

        Parameters
        -----------

        u : MDAnalysis.Universe
            An user-provided MDAnalysis trajectory.

        start : class.int
            frame to start calculation

        stop : class.int
            frame to stop calculation
        
        step : class.int
            number of frames to skip

        Returns
        --------
 
        frames : np.array
            np.array with the selected frames
            
        """

        import numpy as np
        import MDAnalysis as mda

        return np.array(list(map(lambda ts: ts.frame, u.trajectory[start:stop:step])))


def distance_matrix_calc(u, group1_positions, group2_positions):
        """Function to calculate the distane between two atom groups
 
        Parameters
        -----------

        u : MDAnalysis.Universe
            An user-provided MDAnalysis trajectory.

        group1_positions : np.array
            Array with the positions of the first atom group
            
        group2_positions : np.array
            Array with the positions of the second atom group

        Returns
        --------
 
        distance_matrix: np.array
            Matrix with the pair-wise distance between every atoms in group1 and group2

        """

        import MDAnalysis as mda
        from MDAnalysis.lib import distances

        distance_matrix = mda.lib.distances.distance_array(group1_positions, group2_positions, 
                                                           box=u.dimensions, backend="OpenMP")
    

        return distance_matrix


def find_values(array, pairs, len1, len2):
        """Function to find the minors of a matrix
 

        Parameters
        -----------

        array: np.array
            Matrix to perform the calculations on 

        pairs : class.list
            Column and row of specific value to be found in matrix
            
        len1 : class.int
            Length of group 1

        len2 : class.int
            Length of group 2

        Returns
        --------
 
        mat: np.array
            Minor of array

        """
        c = array[pairs[0]:pairs[0]+len1]
        mat = [row[pairs[1]:pairs[1]+len2] for row in c] #minors

        return mat




def get_index_pairs(len1, step1, len2, step2):
        """Function to get atom pairs
 

        Parameters
        -----------

        len1: class.int
            Length of group 1

        step1 : class.int
            Number of atoms of molecules belonging to group 1
            

        len2 : class.int
            Length of group 2

        step2 : class.int
            Number of atom of molecules belonging to group 2


        Returns
        --------
 
        pairs: class.list
            List with atom indices pairs between the two groups

        """

        import numpy as np
        import itertools

        indices1 = np.arange(0, len1, step1)
        indices2 = np.arange(0, len2, step2)
        pairs = filter(lambda x: x[0] != x[1], itertools.product(indices1, indices2))

        return list(pairs)



def coord_matrix_calc(u, atom_name, polymer_oxygen_pos, water_positions, cut_off):
        """Function to calculate coordination number of the chosen atoms of the polymer and water

        Parameters
        -----------

        atom_name: class.str
            polymer atom names to calculate the hydration for

        polymer_oxygen_pos : np.array
            3D np.array with positions of atom_name
           
        water_positions : np.array
            np.array with the positions of the waters


        cut_off : class.float
            number cut off of the distance between the polymer beads and the water for hydration calculation

        Returns
        --------
 
        sums_coord: class.list
            coordination numbers of atom_name

        """
        import numpy as np
        import MDAnalysis as mda
        from MDAnalysis.lib import distances

        group_size=len(atom_name)
       
        distance_matrix = mda.lib.distances.distance_array(polymer_oxygen_pos, water_positions, box= u.dimensions, backend="OpenMP")
        coord_number = np.sum(distance_matrix < cut_off, axis=1)
        sums_coord = [sum(coord_number[i] for i in range(start, len(coord_number), 12)) for start in range(group_size)]
    
        return sums_coord


def pair(array):
        """Function to divde array into pairs
    
        Parameters
        -----------
    
        array: MDAnalysis.Object
          An user-provided array. 
    
        Returns
        --------
    
        a,b: tuple
          A tuple containing the values of the provided array.

        """

        a, b = array[0], array[1]

        return a, b


def get_atom_name_list(atom_sel):
        """Function to obtain the residues combination based on user-provided atom selection.

        Parameter
        ----------

        atom_sel: class.str
          User Provided name/tag of the molecules inside the box.

        Return
        -------

        final: class.list
          List containing all names of atoms.

        """
        import itertools
        
        b = [[str(item), str(item)] for item in atom_sel]
        c = sum(b, [])

        #creating permutation list
        final=list(set(list(itertools.combinations(c, 2)))) 

        return final


def create_network(list_of_tuples):
        """Function to create a network provided by the list of tuples.

        Parameters
        -----------

        list_of_tuples: class.tuple
           List of all edges 

        Returns
        --------
    
        final: networkx.object
           Networkx.object with connected components

        """
        import networkx as nx

        G=nx.Graph()   
        G.add_edges_from(list_of_tuples)
        final=nx.connected_components(G)

        return final


def get_universe_per_frame(u, start, stop, skip):
        """Function to split an MDAnalysis Universe into
        user-provided timeframes.
                
        Parameters
        -----------

        u : MDAnalysis.Universe
            An user-provided MDAnalysis trajectory.

        start : class.int
            frame to start calculation

        stop : class.int
            frame to stop calculation
        
        step : class.int
            number of frames to skip

        Returns
        --------
 
        frames : class.list
            A list with the selected frames

        """

        import MDAnalysis as mda
        
        u_f=[u for ts in u.trajectory[start:stop:skip]]
        #import mda

       

        #u_f=np.array(list(map(lambda ts: ts.frame, u.trajectory[start:stop:skip])))

        return u_f


def create_binned_space(binned_space):


        """Funcion to create a binned array
 

       Parameters
       -----------

        binned_space: np.array
            array with all the bin values


       Returns
       --------
 
       binned_space_volume: np.array
            array with  the bins

        """

        import numpy as np

        binned_space_volume =  list(map(lambda i: [binned_space[i], binned_space[i + 1]], range(len(binned_space) - 1)))

        return binned_space_volume

