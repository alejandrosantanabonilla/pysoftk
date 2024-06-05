Tutorial on the usage of RSA tool
=================================

This tutorial illsutrates how to use the RSA tool on proteins. This
tutorial is divided into two parts. The first part shows how to edit a
protein structure file so that the RSA tool can be ran. The second part
is the running of the RSA tool.

Before starting any analysis, load the neccesary modules for this class.

.. code:: ipython3
    
    from pysoftk.pol_analysis.tools.utils_mda import MDA_input
    from pysoftk.pol_analysis.tools.utils_tools import *
    from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
    from pysoftk.pol_analysis.ring_ring import RSA
    
    import numpy as np
    import pandas as pd
    import MDAnalysis as mda


.. parsed-literal::

    /home/raquellrdc/.local/lib/python3.10/site-packages/tqdm/auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
      from .autonotebook import tqdm as notebook_tqdm


1. Editing the protein structure file
-------------------------------------

Load the initial trajectory

.. code:: ipython3

    topology='data/protein_rsa_noedit.pdb'
    
    trajectory='data/trajectory_rsa_protein.xtc'
    
    u=mda.Universe(topology, trajectory)

.. code:: ipython3

    u.trajectory[0].time

.. parsed-literal::

    150000.0

Now, we need to make sure that each protein can be identified by a
unique resid. That is to say, all atoms of the same protein must have
the same resid, and this resid must be different to the resid of other
proteins in the system. The following code shows how this can be
achieved:

.. code:: ipython3

    proteins_len=len(u.atoms)
    protein_pos=u.atoms.positions
    proteins=u.select_atoms('protein')
    
    proteins.positions=protein_pos
    
    
    #this is a manual way to select the proteins of a system, you need to adapt this code to your own system 
    protein_1=u.select_atoms('segid A')
    protein_2=u.select_atoms('segid B')
    protein_3=u.select_atoms('segid C')
    protein_4=u.select_atoms('segid D')
    protein_5=u.select_atoms('segid E')
    protein_6=u.select_atoms('segid F')
    protein_7=u.select_atoms('segid G')
    protein_8=u.select_atoms('segid H')
    protein_9=u.select_atoms('segid I')
    protein_10=u.select_atoms('segid J')
    
    
    protein_1_len=len(np.unique(protein_1.resids))
    protein_2_len=len(np.unique(protein_2.resids))
    protein_3_len=len(np.unique(protein_3.resids))
    protein_4_len=len(np.unique(protein_4.resids))
    protein_5_len=len(np.unique(protein_5.resids))
    protein_6_len=len(np.unique(protein_6.resids))
    protein_7_len=len(np.unique(protein_7.resids))
    protein_8_len=len(np.unique(protein_8.resids))
    protein_9_len=len(np.unique(protein_9.resids))
    protein_10_len=len(np.unique(protein_10.resids))
    
    print(protein_1_len)
    
    
    
    
    resids=[]
    
    
    resids= resids + [1]*protein_1_len
    resids= resids + [2]*protein_2_len
    resids= resids + [3]*protein_3_len
    
    resids= resids + [4]*protein_4_len
    resids= resids + [5]*protein_5_len
    resids= resids + [6]*protein_6_len
    
    resids= resids + [7]*protein_7_len
    resids= resids + [8]*protein_8_len
    resids= resids + [9]*protein_9_len
    
    resids= resids + [10]*protein_10_len
    
    
    
    print(len(resids))
    print(len(proteins))
    
    proteins.residues.resids=resids
    
    with mda.Writer("data/trajectory_resids.pdb", proteins.n_atoms) as W:
        
        W.write(proteins)


.. parsed-literal::

    31
    310
    3270


Great, now the new structure file has been created, so that each protein
has a unique resid for all of its atoms. Now we can run the RSA tool on
this new file.

2.Running RSA
-------------

With the new structure file, we can run the RSA tool in the exact same
way we run it on polymers

.. code:: ipython3

    
    topology='data/trajectory_resids.pdb'
    
    trajectory='data/trajectory_rsa_protein.xtc'

.. code:: ipython3

    #name output file
    results='data/rsa_prot_tutorial.parquet'
    
    #angle cutoff - angle range (val < ang_c or val> 180-ang_c). 
    ang_c=30
    
    
    #distance cutoff - distance between two rings to be considered stacked
    dist_c=5
    


.. code:: ipython3

    rsa=RSA(topology, trajectory).stacking_analysis(dist_c, ang_c, 0, 20, 2, results)


.. parsed-literal::

    Ring Stacking analysis has started


.. parsed-literal::

    100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 45/45 [00:00<00:00, 548.29it/s]
    Detecting atoms in Rings: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:24<00:00,  1.08it/s]
    Separating Rings: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:00<00:00, 11647.11it/s]


.. parsed-literal::

    
    
    Preparing DataFrame to store the results


.. parsed-literal::

      0%|                                                                                                                                                                                      | 0/10 [00:00<?, ?it/s]
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:08,  2.91it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:00<00:04,  5.06it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:00<00:04,  4.55it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:01<00:04,  4.97it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:01<00:04,  3.85it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:05,  3.51it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:05,  3.00it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:02<00:05,  2.81it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:03<00:04,  3.30it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:03<00:04,  2.99it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:04<00:04,  2.77it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:04<00:03,  2.79it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:03,  2.59it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:05<00:03,  2.58it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:05<00:03,  2.55it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:06<00:02,  3.14it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:06<00:01,  3.41it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:06<00:01,  2.89it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:07<00:01,  3.00it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:07<00:01,  2.76it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:07<00:00,  2.79it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:08<00:00,  2.63it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:08<00:00,  3.02it/s][A
     10%|█████████████████▍                                                                                                                                                            | 1/10 [00:08<01:17,  8.61s/it]

.. parsed-literal::

    ([(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([2173, 2174, 2177, 2179, 2182]), array([2411, 2413, 2403, 2404, 2406, 2408]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831])), (array([2500, 2501, 2504, 2506, 2509]), array([2706, 2708, 2711, 2702, 2703]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:04,  5.93it/s][A
    Computing stacking distances:   8%|███████████                                                                                                                                     | 2/26 [00:00<00:05,  4.11it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:00<00:07,  2.93it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:08,  2.58it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:05,  3.50it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:01<00:05,  3.60it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:01<00:03,  5.68it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:04,  4.23it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:02<00:03,  4.99it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:02<00:03,  4.51it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:03<00:04,  3.45it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:03<00:03,  4.07it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:03<00:03,  3.34it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:04<00:03,  3.03it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:04<00:03,  2.84it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:04<00:03,  2.64it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:05<00:03,  2.52it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:05<00:02,  2.45it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:06<00:02,  2.91it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:06<00:01,  2.81it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:06<00:01,  3.50it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:06<00:00,  3.51it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:07<00:00,  2.93it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:07<00:00,  2.64it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:08<00:00,  3.21it/s][A
     20%|██████████████████████████████████▊                                                                                                                                           | 2/10 [00:16<01:06,  8.31s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([193, 195, 198, 200, 190, 191]), array([871, 874, 865, 866, 869]))], [(array([211, 212, 215, 217, 220]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))], [(array([1749, 1750, 1752, 1754, 1757, 1759]), array([3133, 3134, 3136, 3138, 3141, 3143]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831]))], [(array([2809, 2811, 2814, 2816, 2806, 2807]), array([3065, 3067, 3057, 3058, 3060, 3062]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:11,  2.24it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:00<00:04,  4.79it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:06,  3.60it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:06,  3.01it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:01<00:05,  3.63it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:06,  3.02it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:06,  2.77it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:04,  3.47it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:02<00:04,  3.70it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:03<00:04,  3.02it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:03<00:04,  3.08it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:04<00:04,  2.77it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:04<00:04,  2.47it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:03,  2.87it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:05<00:03,  2.65it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:06<00:03,  2.41it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:06<00:02,  2.35it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:06<00:02,  2.87it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:07<00:01,  3.09it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:07<00:01,  3.04it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:07<00:01,  2.81it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:08<00:00,  2.80it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:08<00:00,  2.57it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:09<00:00,  2.86it/s][A
     30%|████████████████████████████████████████████████████▏                                                                                                                         | 3/10 [00:25<01:00,  8.67s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([847, 849, 852, 854, 844, 845]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1067, 1068, 1071, 1073, 1076]), array([1852, 1855, 1846, 1847, 1850]))], [(array([2173, 2174, 2177, 2179, 2182]), array([2411, 2413, 2403, 2404, 2406, 2408]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831])), (array([2500, 2501, 2504, 2506, 2509]), array([2706, 2708, 2711, 2702, 2703]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:10,  2.46it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:00<00:05,  3.86it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:06,  3.27it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:06,  3.49it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:01<00:05,  3.50it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:06,  3.15it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:05,  3.07it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:05,  3.14it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:03<00:05,  2.90it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:03<00:05,  2.95it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:03<00:05,  2.70it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:04<00:04,  2.71it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:04<00:04,  2.70it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:05<00:04,  2.56it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:04,  2.39it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:06<00:03,  2.32it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:06<00:03,  2.48it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:06<00:02,  2.39it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:07<00:02,  2.36it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:07<00:01,  2.85it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:07<00:01,  3.41it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:08<00:01,  2.96it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:08<00:00,  2.83it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:08<00:00,  2.89it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:09<00:00,  2.82it/s][A
     40%|█████████████████████████████████████████████████████████████████████▌                                                                                                        | 4/10 [00:35<00:53,  8.89s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))], [(array([1192, 1193, 1196, 1198, 1201]), array([3160, 3163, 3154, 3155, 3158]))], [(array([2173, 2174, 2177, 2179, 2182]), array([2411, 2413, 2403, 2404, 2406, 2408]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:11,  2.15it/s][A
    Computing stacking distances:   8%|███████████                                                                                                                                     | 2/26 [00:00<00:10,  2.23it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:01<00:09,  2.48it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:07,  2.75it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:02<00:08,  2.50it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:02<00:08,  2.49it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:07,  2.43it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:03<00:06,  2.63it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:03<00:05,  3.39it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:03<00:05,  2.85it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:04<00:05,  2.99it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:04<00:04,  3.30it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:04<00:04,  2.91it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:05<00:04,  2.68it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:02,  4.37it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:05<00:02,  3.97it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:05<00:01,  4.52it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:06<00:01,  4.23it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:06<00:01,  3.60it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:07<00:01,  3.11it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:07<00:00,  3.46it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:07<00:00,  2.93it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:08<00:00,  2.73it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:08<00:00,  3.13it/s][A
     50%|███████████████████████████████████████████████████████████████████████████████████████                                                                                       | 5/10 [00:43<00:43,  8.68s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([193, 195, 198, 200, 190, 191]), array([871, 874, 865, 866, 869]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([847, 849, 852, 854, 844, 845]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1067, 1068, 1071, 1073, 1076]), array([1852, 1855, 1846, 1847, 1850])), (array([1192, 1193, 1196, 1198, 1201]), array([1725, 1727, 1730, 1721, 1722]))], [(array([2173, 2174, 2177, 2179, 2182]), array([2411, 2413, 2403, 2404, 2406, 2408]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:11,  2.16it/s][A
    Computing stacking distances:   8%|███████████                                                                                                                                     | 2/26 [00:00<00:09,  2.54it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:01<00:09,  2.34it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:08,  2.69it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:07,  2.94it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:04,  4.73it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:04,  3.75it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:03,  4.44it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:02<00:04,  3.62it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:03<00:02,  5.47it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:03<00:02,  4.65it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:03<00:03,  3.60it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:04<00:03,  3.28it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:04<00:03,  2.89it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:04<00:02,  3.17it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:05<00:02,  3.56it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:05<00:02,  3.04it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:06<00:02,  2.81it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:06<00:01,  3.29it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:06<00:01,  2.97it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:06<00:01,  2.90it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:07<00:00,  3.00it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:07<00:00,  3.34it/s][A
     60%|████████████████████████████████████████████████████████████████████████████████████████████████████████▍                                                                     | 6/10 [00:51<00:33,  8.38s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([211, 212, 215, 217, 220]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([1067, 1068, 1071, 1073, 1076]), array([1852, 1855, 1846, 1847, 1850])), (array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:10,  2.39it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:00<00:07,  3.20it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:07,  2.81it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:01<00:06,  3.29it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:06,  2.98it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:05,  3.02it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:05,  3.27it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:03<00:05,  3.09it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:03<00:05,  2.75it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:03<00:04,  3.06it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:04<00:04,  2.78it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:04<00:03,  3.25it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:03,  2.98it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:05<00:03,  2.80it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:06<00:02,  2.98it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:06<00:02,  3.03it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:06<00:02,  2.71it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:07<00:01,  2.53it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:07<00:01,  2.42it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:08<00:01,  2.27it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:08<00:00,  2.19it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:09<00:00,  2.13it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:09<00:00,  2.69it/s][A
     70%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▊                                                    | 7/10 [01:00<00:26,  8.81s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([847, 849, 852, 854, 844, 845]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:09,  2.55it/s][A
    Computing stacking distances:   8%|███████████                                                                                                                                     | 2/26 [00:00<00:07,  3.03it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:01<00:08,  2.84it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:08,  2.55it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:08,  2.43it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:02<00:07,  2.53it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:06,  3.05it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:05,  3.21it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:02<00:04,  3.58it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:03<00:05,  3.03it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:03<00:05,  2.75it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:04<00:05,  2.60it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:04<00:04,  2.61it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:04<00:03,  3.15it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:05<00:03,  2.81it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:03,  2.77it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:06<00:03,  2.68it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:06<00:03,  2.53it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:06<00:02,  2.42it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:07<00:02,  2.36it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:07<00:01,  2.99it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:08<00:01,  2.69it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:08<00:01,  2.63it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:08<00:00,  3.26it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:08<00:00,  3.18it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:09<00:00,  2.80it/s][A
     80%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▏                                  | 8/10 [01:10<00:17,  8.96s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([847, 849, 852, 854, 844, 845]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1067, 1068, 1071, 1073, 1076]), array([1852, 1855, 1846, 1847, 1850])), (array([1174, 1176, 1179, 1181, 1171, 1172]), array([1852, 1855, 1846, 1847, 1850])), (array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))], [(array([1192, 1193, 1196, 1198, 1201]), array([3160, 3163, 3154, 3155, 3158]))], [(array([2173, 2174, 2177, 2179, 2182]), array([2411, 2413, 2403, 2404, 2406, 2408]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:08,  2.93it/s][A
    Computing stacking distances:   8%|███████████                                                                                                                                     | 2/26 [00:00<00:06,  3.64it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:01<00:07,  2.88it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:08,  2.64it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:07,  2.67it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:02<00:07,  2.54it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:07,  2.48it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:03<00:07,  2.44it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:03<00:06,  2.44it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:03<00:06,  2.34it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:04<00:06,  2.30it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:04<00:05,  2.77it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:05<00:04,  2.64it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:05<00:04,  2.59it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:05<00:03,  3.17it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:03,  3.24it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:06<00:03,  2.92it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:06<00:02,  2.97it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:07<00:02,  2.65it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:07<00:02,  2.45it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:08<00:02,  2.34it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:08<00:01,  2.48it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:08<00:00,  2.92it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:09<00:00,  2.79it/s][A
     90%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                 | 9/10 [01:19<00:09,  9.07s/it]

.. parsed-literal::

    ([(array([86, 87, 90, 92, 95]), array([544, 547, 538, 539, 542]))], [(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([847, 849, 852, 854, 844, 845]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1095, 1096, 1098, 1100, 1103, 1105]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1067, 1068, 1071, 1073, 1076]), array([1852, 1855, 1846, 1847, 1850])), (array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))], [(array([1192, 1193, 1196, 1198, 1201]), array([3160, 3163, 3154, 3155, 3158]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831]))])


.. parsed-literal::

    
    Computing stacking distances:   0%|                                                                                                                                                        | 0/26 [00:00<?, ?it/s][A
    Computing stacking distances:   4%|█████▌                                                                                                                                          | 1/26 [00:00<00:11,  2.18it/s][A
    Computing stacking distances:   8%|███████████                                                                                                                                     | 2/26 [00:00<00:11,  2.18it/s][A
    Computing stacking distances:  12%|████████████████▌                                                                                                                               | 3/26 [00:01<00:10,  2.24it/s][A
    Computing stacking distances:  15%|██████████████████████▏                                                                                                                         | 4/26 [00:01<00:07,  2.95it/s][A
    Computing stacking distances:  19%|███████████████████████████▋                                                                                                                    | 5/26 [00:01<00:06,  3.34it/s][A
    Computing stacking distances:  23%|█████████████████████████████████▏                                                                                                              | 6/26 [00:02<00:05,  3.35it/s][A
    Computing stacking distances:  27%|██████████████████████████████████████▊                                                                                                         | 7/26 [00:02<00:05,  3.31it/s][A
    Computing stacking distances:  31%|████████████████████████████████████████████▎                                                                                                   | 8/26 [00:02<00:06,  2.82it/s][A
    Computing stacking distances:  35%|█████████████████████████████████████████████████▊                                                                                              | 9/26 [00:03<00:06,  2.63it/s][A
    Computing stacking distances:  38%|███████████████████████████████████████████████████████                                                                                        | 10/26 [00:03<00:05,  2.90it/s][A
    Computing stacking distances:  42%|████████████████████████████████████████████████████████████▌                                                                                  | 11/26 [00:03<00:05,  2.66it/s][A
    Computing stacking distances:  46%|██████████████████████████████████████████████████████████████████                                                                             | 12/26 [00:04<00:05,  2.62it/s][A
    Computing stacking distances:  50%|███████████████████████████████████████████████████████████████████████▌                                                                       | 13/26 [00:04<00:05,  2.45it/s][A
    Computing stacking distances:  54%|█████████████████████████████████████████████████████████████████████████████                                                                  | 14/26 [00:05<00:05,  2.28it/s][A
    Computing stacking distances:  58%|██████████████████████████████████████████████████████████████████████████████████▍                                                            | 15/26 [00:05<00:03,  2.77it/s][A
    Computing stacking distances:  62%|████████████████████████████████████████████████████████████████████████████████████████                                                       | 16/26 [00:05<00:03,  2.93it/s][A
    Computing stacking distances:  65%|█████████████████████████████████████████████████████████████████████████████████████████████▌                                                 | 17/26 [00:06<00:02,  3.26it/s][A
    Computing stacking distances:  69%|███████████████████████████████████████████████████████████████████████████████████████████████████                                            | 18/26 [00:06<00:02,  2.94it/s][A
    Computing stacking distances:  73%|████████████████████████████████████████████████████████████████████████████████████████████████████████▌                                      | 19/26 [00:06<00:02,  3.11it/s][A
    Computing stacking distances:  77%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████                                 | 20/26 [00:07<00:01,  3.06it/s][A
    Computing stacking distances:  81%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                           | 21/26 [00:07<00:01,  2.69it/s][A
    Computing stacking distances:  85%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████                      | 22/26 [00:07<00:01,  3.13it/s][A
    Computing stacking distances:  88%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌                | 23/26 [00:08<00:01,  2.79it/s][A
    Computing stacking distances:  92%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████           | 24/26 [00:08<00:00,  2.56it/s][A
    Computing stacking distances:  96%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▌     | 25/26 [00:09<00:00,  2.48it/s][A
    Computing stacking distances: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 26/26 [00:09<00:00,  2.74it/s][A
    100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 10/10 [01:28<00:00,  8.89s/it]

.. parsed-literal::

    ([(array([740, 741, 744, 746, 749]), array([1171, 1172, 1174, 1176, 1179, 1181]))], [(array([847, 849, 852, 854, 844, 845]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1095, 1096, 1098, 1100, 1103, 1105]), array([1498, 1499, 1501, 1503, 1506, 1508]))], [(array([1192, 1193, 1196, 1198, 1201]), array([1852, 1855, 1846, 1847, 1850]))], [(array([1192, 1193, 1196, 1198, 1201]), array([3160, 3163, 3154, 3155, 3158]))], [(array([1749, 1750, 1752, 1754, 1757, 1759]), array([3133, 3134, 3136, 3138, 3141, 3143]))], [(array([2173, 2174, 2177, 2179, 2182]), array([2411, 2413, 2403, 2404, 2406, 2408]))], [(array([2482, 2484, 2487, 2489, 2479, 2480]), array([2833, 2836, 2827, 2828, 2831]))])
    Information succesfully stored in data/rsa_prot_tutorial.parquet
    Stacking analysis has succesfully finished!
    Function stacking_analysis Took 113.2266 seconds


.. code:: ipython3

    
    df_results = 'data/rsa_prot_tutorial.parquet'
    df = pd.read_parquet(df_results)
    
    df




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>atom_index</th>
          <th>pol_resid</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>[[[740, 741, 744, 746, 749], [1171, 1172, 1174...</td>
          <td>[[[3, 4], [7, 8], [8, 9]]]</td>
        </tr>
        <tr>
          <th>1</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [1, 3], [1, 4], [3, 4], [4, 6], [6, ...</td>
        </tr>
        <tr>
          <th>2</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [3, 4], [3, 5], [4, 6], [7, 8], [8, ...</td>
        </tr>
        <tr>
          <th>3</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [3, 4], [4, 6], [4, 10], [7, 8]]]</td>
        </tr>
        <tr>
          <th>4</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [1, 3], [3, 4], [3, 5], [4, 6], [7, ...</td>
        </tr>
        <tr>
          <th>5</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [1, 4], [3, 4], [4, 6]]]</td>
        </tr>
        <tr>
          <th>6</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [3, 4], [3, 5], [4, 6], [8, 9]]]</td>
        </tr>
        <tr>
          <th>7</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [3, 4], [3, 5], [4, 6], [4, 10], [7,...</td>
        </tr>
        <tr>
          <th>8</th>
          <td>[[[86, 87, 90, 92, 95], [544, 547, 538, 539, 5...</td>
          <td>[[[1, 2], [3, 4], [3, 5], [4, 5], [4, 6], [4, ...</td>
        </tr>
        <tr>
          <th>9</th>
          <td>[[[740, 741, 744, 746, 749], [1171, 1172, 1174...</td>
          <td>[[[3, 4], [3, 5], [4, 5], [4, 6], [4, 10], [6,...</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    df.iloc[0][0]




.. parsed-literal::

    array([array([array([740, 741, 744, 746, 749]),
                  array([1171, 1172, 1174, 1176, 1179, 1181]),
                  array([2173, 2174, 2177, 2179, 2182]),
                  array([2411, 2413, 2403, 2404, 2406, 2408]),
                  array([2482, 2484, 2487, 2489, 2479, 2480]),
                  array([2833, 2836, 2827, 2828, 2831]),
                  array([2500, 2501, 2504, 2506, 2509]),
                  array([2706, 2708, 2711, 2702, 2703])], dtype=object)],
          dtype=object)



.. code:: ipython3

    sev_ring=RSA(topology, trajectory).find_several_rings_stacked(df_results)

.. code:: ipython3

    #print the resids of the networ of polymers connected by their stacked rings
    print(sev_ring)


.. parsed-literal::

    [[{3, 4}, {8, 9, 7}], [{1, 2, 3, 4, 6, 8, 9, 10}], [{1, 2}, {3, 4, 5, 6}, {8, 9, 7}], [{1, 2}, {10, 3, 4, 6}, {8, 7}], [{1, 2, 3, 4, 5, 6}, {8, 9, 7}], [{1, 2, 3, 4, 6}], [{1, 2}, {3, 4, 5, 6}, {8, 9}], [{1, 2}, {3, 4, 5, 6, 10}, {8, 7}], [{1, 2}, {3, 4, 5, 6, 10}, {8, 9}], [{3, 4, 5, 6, 10}, {8, 9, 7}]]


