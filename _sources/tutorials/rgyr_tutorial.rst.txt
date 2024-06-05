Tutorial on the usage of rgyr_micelle
=====================================

This tutorial will show how to use the rgyr_micelle method to compute
the radius of gyration of a micelle that is broken across the periodic
boundary conditions (PBC) at some time steps. It is worth noticing that
this class can be used in a similar manner for aggregates that are
connected across the PBC

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
    from rgyr_micelle import rgyr
    #from pysoftk.pol_analysis.rgyr_micelle import rgyr
    
    import numpy as np
    import pandas as pd


.. parsed-literal::

    /home/raquellrdc/Desktop/PhD/mda_umap/lib/python3.10/site-packages/tqdm/auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
      from .autonotebook import tqdm as notebook_tqdm


1. First load the trajectory, we are going to use the cyclic topology
   used in previous tutorials.

.. code:: ipython3

    topology='data/short_movie_cyclic.tpr'
    trajectory='data/short_movie_cyclic.xtc'


2. Using SCP, we need to cluster the polymers of the simulations, to
   obtain the resids of the aggregate that we want to perform the
   analysis on. Following the same steps as in the SCP_tutorial for the
   cyclic topology.

.. code:: ipython3

    results_name='data/pictures_tutorial/cyclic_scp_result'
    
    atom_names=['C02T']
    
    cluster_cutoff = 12
    
    start=0
    stop=10001
    step=1
    
    
    #running SCP
    c = SCP(topology, trajectory).spatial_clustering_run(start, stop, step, atom_names, cluster_cutoff, results_name)


.. parsed-literal::

    100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00, 46.12it/s]

.. parsed-literal::

    The file data/pictures_tutorial/cyclic_scp_result.parquet has been successfully created.
    Function spatial_clustering_run Took 0.1947 seconds


3. Now that the polymers are clustered, we are going to select the
   largest aggregate to perform the analysis on. We can use
   make_micelle_whole to obtain the resids of this larger aggregates for
   all time steps. Note that you can select whichever aggregate by
   inspecting the clustering pandas dataframe result.

.. code:: ipython3

    
    #the result from SCP
    resids_total='data/pictures_tutorial/cyclic_scp_result.parquet'
    
    #obtaining the largest aggregate
    largest_micelle_resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(resids_total)

4. Let’s obtain the whole positions of the largest micelle with
   micelle_whole. This is done in the same way as in the micelle_whole
   tutorial.

.. code:: ipython3

    #select the resname of the polymers
    resname=['LIG'] 
    
    #run to obtain the whole positions
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(resname, largest_micelle_resids, start, stop, step)


.. parsed-literal::

      0%|                                                                                                                                                    | 0/3 [00:00<?, ?it/s]/home/raquellrdc/Desktop/PhD/pysoftk/alejandro_newest_releast_check/pysoftk_analysis_code/test_final/make_micelle_whole.py:347: FutureWarning: arrays to stack must be passed as a "sequence" type such as list or tuple. Support for non-sequence iterables such as generators is deprecated as of NumPy 1.16 and will raise an error in the future.
      atom_positions_over_trajectory = list(tqdm(map(self.make_cluster_whole, frames, resname, cluster_resids_f[0],
    100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00,  7.22it/s]

.. parsed-literal::

    Elapsed time for matrix calculation: 0.4943 seconds


5. Now, we have all the neccesary inputs to calculate the radius of
   gyration of the micelle. All that we need is the resids of the
   polymers that belong to the micelle that we want to perform the
   analysis of, their whole positions and the frames on which we want to
   run the analysis on. Note that these frames need to be the same as in
   the make_micelle_whole calculation.

.. code:: ipython3

    rgyr_micelle_whole = rgyr(topology, trajectory).running_rgyr(largest_micelle_resids, atom_pos, start, stop, step)



.. parsed-literal::

    100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00, 1451.15it/s]


And this is the radius of gyration taken the PBC properly into account!
