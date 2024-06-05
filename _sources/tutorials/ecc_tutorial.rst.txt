Tutorial on the usage of ecc_micelle
====================================

This tutorial shows how to use the ecc_micelle class to compute the
eccentricity of a micelle that is broken across the periodic boundary
conditions (PBC) at some time steps. Please note that this class can be
used in the same way for aggregates that are not broken across the PBC.

Before starting any analysis, load the neccesary modules for this class.

.. code:: ipython3

    from  utils_mda import MDA_input
    #from pysoftk.pol_analysis.tools.utils_mda import MDA_input
    from utils_tools import *
    #from pysoftk.pol_analysis.tools.utils_tools import *
    from clustering import SCP
    #from pysoftk.pol_analysis.clustering import SCP
    from make_micelle_whole import micelle_whole
    #from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
    from ecc_micelle import ecc
    #from pysoftk.pol_analysis.ecc_micelle import ecc
    
    import numpy as np
    import pandas as pd

1. First load the trajectory, we are going to use the cyclic topology
   used in previous tutorials.

.. code:: ipython3

    topology='data/short_movie_cyclic.tpr'
    trajectory='data/short_movie_cyclic.xtc'


2. Import the clustering data from SCP function

.. code:: ipython3

    resids_total='data/pictures_tutorial/cyclic_scp_result.parquet'

3. Obtain the largest micelle from the clustering pandas dataframe

.. code:: ipython3

    largest_micelle_resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(resids_total)

4. Now, obtain the positions of the whole micelle

.. code:: ipython3

    #select the resname of the polymers
    resname=['LIG']
    
    #run to obtain the whole positions
    start=0
    stop=10001
    step=1
    
    #run micelle_whole
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(resname, largest_micelle_resids, start, stop, step)


.. parsed-literal::

      0%|                                                                                                                                                    | 0/3 [00:00<?, ?it/s]/home/raquellrdc/Desktop/PhD/pysoftk/alejandro_newest_releast_check/pysoftk_analysis_code/test_final/make_micelle_whole.py:347: FutureWarning: arrays to stack must be passed as a "sequence" type such as list or tuple. Support for non-sequence iterables such as generators is deprecated as of NumPy 1.16 and will raise an error in the future.
      atom_positions_over_trajectory = list(tqdm(map(self.make_cluster_whole, frames, resname, cluster_resids_f[0],
    100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00,  6.86it/s]

.. parsed-literal::

    Elapsed time for matrix calculation: 0.5386 seconds

5. We have all the neccesary inputs to run the eccentricity calculation

.. code:: ipython3

    ecc_micelle_whole = ecc(topology, trajectory).running_ecc(largest_micelle_resids, atom_pos, start, stop, step)
    ecc_micelle_whole


.. parsed-literal::

    100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00, 1241.41it/s]




.. parsed-literal::

    [0.31216114645263193, 0.09429518042597795, 0.0661655070142777]



That is it! ecc_micelle_whole contains the eccentricity of the micelle
for the selected timesteps
