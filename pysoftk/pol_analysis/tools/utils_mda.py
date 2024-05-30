import MDAnalysis as mda


class MDA_input(object):
    """ Class to load in memory a trajectory using MDA.

    Attributes
    ----------

    xtc : str
       Path and Name of the file with extension xtc.

    tpr : str
       Path and Name of the file with extension tpr.

    """

    def __init__(self, tpr_file, xtc_file):
        """
        
        Parameters
        -----------

        xtc : str
          Path and Name of the file with extension xtc.

        tpr : str
          Path and Name of the file with extension tpr.
        
        """  
        self.tpr_file = tpr_file
        self.xtc_file = xtc_file 

    def get_mda_universe(self):
        """Function to create a MDAnalysis object loading 
           the trajectory.

        Parameters
        -----------

        None
        

        Returns
        --------

        u: MDAnalysis.object
          MDAnalysis object containing relevant information.

        """

        u=mda.Universe(self.tpr_file, self.xtc_file)

        return u
    

   