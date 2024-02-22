import os
import time
from functools import wraps

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.tools.utils_tools import *

class ecc(MDA_input):
    """A class used to compute the eccentricity of micelles 
    """
    
    def __init__(self, tpr_file, xtc_file):
      """  Instantiating the MDA_input class as
          super function.
      """
      
      super().__init__(tpr_file, xtc_file)


    def ecc_calc(self, atom_sel_micelle, micelle_positions):
        """Function to calculate the eccentricity of the chosen atoms

        Parameters
        -----------

        atom_sel_micelle : MDAnalysis.Object
            MDAnalysis.atom_selection function
            

        micelle_positions : np.array
            np.array with the positions of the micelle atoms

        Returns
        --------
 
        ecc: class.float
            Eccentricity of the selected atoms

        """

        import MDAnalysis as mda
        import numpy as np

        micelle = atom_sel_micelle
        micelle.positions = micelle_positions[:][1]


        moment_of_inertia = (micelle.moment_of_inertia()).diagonal()  
        ecc= 1-np.min(moment_of_inertia)/np.mean(moment_of_inertia)   

        return ecc
    

    def running_ecc(self, atom_names, micelle_positions, start, stop, skip):
        """Function to calculate the eccentricity of the chosen atoms over time
 
        Parameters
        -----------

        atom_sel_type : class.str
            Type of atom selection

        atom_names : class.list
            List with the names of atoms
            
        micelle_positions : np.array
            np.array with the positions of the micelle atoms

        Start : class.int
            Starting frame of the trajectory

        Stop : class.int
            Ending rame of the trajectory

        Skip : class.int
            Skipping every that many framesof the trajectory

        Returns
        --------
 
        ecc: np.array
            Eccentricity of the selected atoms over time

        """

        import MDAnalysis as mda
        from tqdm.auto import tqdm
        

        u=super().get_mda_universe()
     
        us = [u]*len(u.trajectory[start:stop:skip])        
        atom_sel_micelle = map(selecting_atoms, us, atom_names)
        
        ecc_f = list(tqdm(map(self.ecc_calc, atom_sel_micelle, micelle_positions), total=len(us)))

        return ecc_f
