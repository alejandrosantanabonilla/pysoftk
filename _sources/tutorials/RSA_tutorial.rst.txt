Tutorial on the usage of RSA tool
=================================

This tutorial illsutrates how to use the RSA tool.

Before starting any analysis, load the neccesary modules for this class.

.. code:: ipython3

    
    from pysoftk.pol_analysis.tools.utils_mda import MDA_input
    from pysoftk.pol_analysis.tools.utils_tools import *
    from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
    from pysoftk.pol_analysis.ring_ring import RSA
    
    import numpy as np
    import pandas as pd


.. parsed-literal::

    /home/raquellrdc/.local/lib/python3.10/site-packages/tqdm/auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
      from .autonotebook import tqdm as notebook_tqdm


1. Select your trajectory files, it is recommended to use a tpr file for
   the topology and xtc file for the trajectory. Note that any
   MDAnalysis supported file can be used here.

.. code:: ipython3

    topology='data/f8bt_slab_quench.tpr'
    trajectory='data/1_frame_traj.xtc'

The simulation where we are going to perform the analysis is on this
very big system. Since the system is very big, we will only perform the
analysis on one frame.

.. figure:: images/ring_system.png
   :alt: Image Alt Text

Which is a simulation box filled with this polymer:

.. figure:: images/ring_1_polymer.png
   :alt: Image Alt Text

This class requires minimal user input. Only the angle, distance cutoff,
start, stop and step frames are needed, as well as the name of the
output file.

.. code:: ipython3

    #name output file
    results='data/rsa.parquet'
    
    #angle cutoff - angle range (val < ang_c or val> 180-ang_c). 
    ang_c=30
    
    #distance cutoff - distance between two rings to be considered stacked
    dist_c=5
    
    #start frame
    start=0
    
    #stop frame
    stop=1
    
    #step frame
    step=1
    
Now, let’s run the RSA stacking analysis!-This run will take a bit of
time, it is mainly because our system is very large!

.. code:: ipython3

    rsa=RSA(topology, trajectory).stacking_analysis(dist_c, ang_c, start, stop, step, results)


.. parsed-literal::

    Ring Stacking analysis has started


.. parsed-literal::

    
      0%|                                                                                                                                                                                  | 0/300700 [00:00<?, ?it/s][A
      0%|                                                                                                                                                                        | 15/300700 [00:00<34:59, 143.20it/s][A
      0%|                                                                                                                                                                        | 30/300700 [00:00<35:02, 143.00it/s][A
      0%|                                                                                                                                                                        | 45/300700 [00:00<35:03, 142.94it/s][A

      .
      .
      .
      
      7%|███████████▎                                                                                                                                                         | 20597/300700 [02:15<34:43, 134.44it/s][A
      7%|███████████▎                                                                                                                                                         | 20612/300700 [02:15<33:46, 138.20it/s][A

.. code:: ipython3

    
    df_results = 'data/rsa.parquet'
    df = pd.read_parquet(df_results)
    print(df)

You will see that a lot of pdb files have been printed in the directory.
These are just some results from the output where you can see visually
the ring stacking! 

.. image:: images/ring_stacking_snapshot_1.png
.. image:: images/ring_stacking_snapshot_2.png
.. image:: images/ring_stacking_snapshot_3.png

We can make use of another function in the RSA class that output the
network of polymers that have their rings stacked. All we need is the
pandas dataframe outputed by the ring stacking calculation

.. code:: ipython3

    sev_ring=RSA(topology, trajectory).find_several_rings_stacked(df_results)

.. code:: ipython3

    #print the resids of the network of polymers connected by their stacked rings
    print(sev_ring)

.. code:: ipython3

    #print the resids of the network of polymers connected by their stacked rings
    print(sev_ring)


.. parsed-literal::

    [{1, 322, 642, 262, 620, 216, 239, 212, 276, 182, 728, 348, 223}, {480, 2, 68, 10, 20, 88, 30, 607}, {96, 4, 558}, {131, 5, 393, 139, 396, 269, 143, 16, 15, 146, 661, 22, 150, 25, 539, 29, 158, 159, 37, 550, 293, 424, 165, 170, 299, 172, 171, 46, 426, 186, 187, 188, 444, 61, 576, 449, 326, 71, 327, 460, 79, 80, 81, 720, 84, 215, 90, 219, 221, 222, 352, 482, 226, 354, 363, 110, 366, 751, 113, 116, 244, 118, 119, 628, 381}, {321, 69, 7, 41, 530, 213, 379}, {514, 9, 457, 398, 441, 414}, {672, 194, 12, 677}, {42, 53, 13}, {504, 14, 431}, {161, 21}, {26, 372}, {207, 31}, {32, 98, 707, 740, 712, 747, 716, 589, 752, 369, 723, 343, 761, 765}, {33, 405, 551}, {34, 534}, {291, 35, 773, 36, 238, 402, 179, 283}, {74, 38}, {533, 39}, {40, 629, 166}, {481, 47}, {353, 134, 136, 776, 75, 50, 274, 277, 311, 760, 58, 315, 157}, {130, 518, 715, 51, 632, 731, 351}, {388, 517, 390, 526, 144, 531, 403, 419, 547, 421, 552, 554, 556, 175, 561, 690, 52, 438, 184, 447, 477, 99, 107, 109, 510}, {120, 185, 557, 54}, {501, 55}, {259, 680, 496, 497, 85, 56, 665, 639}, {57, 82, 173}, {227, 389, 618, 314, 399, 623, 273, 466, 471, 250, 667, 60, 604}, {736, 705, 385, 706, 775, 296, 137, 555, 108, 652, 686, 335, 588, 370, 63}, {64, 432}, {65, 247, 367}, {418, 67, 70, 422, 622, 465}, {72, 493, 83, 440, 89}, {73, 140}, {76, 365, 94}, {261, 519, 104, 77, 270, 529, 564, 411, 574}, {257, 129, 453, 429, 78}, {91, 443}, {768, 169, 92, 595, 700}, {97, 196}, {192, 577, 100, 292}, {640, 512, 645, 774, 391, 670, 439, 704, 452, 583, 458, 208, 596, 597, 724, 612, 106, 763, 636}, {486, 154, 168, 553, 117, 122, 575}, {123, 340, 373, 318}, {342, 127}, {128, 521, 155}, {258, 132, 374, 303}, {320, 325, 133, 744, 302, 469, 633}, {537, 339, 741, 135}, {344, 147, 316, 333}, {153, 713}, {520, 156, 408}, {160, 545}, {163, 317}, {456, 167}, {177, 682, 689, 599}, {416, 407, 549, 178, 246, 572, 473, 508, 190}, {193, 298, 543}, {201, 195}, {307, 371, 197}, {225, 611, 198, 205, 206, 657, 658, 217}, {256, 464, 199}, {200, 264, 210, 721}, {660, 664, 666, 290, 678, 436, 202, 334, 337, 211, 468, 341, 729, 606, 479, 230, 242, 502, 249, 638}, {234, 203}, {204, 294}, {392, 214}, {224, 267, 286}, {323, 265, 361, 235, 309}, {364, 245}, {248, 251, 693, 483}, {252, 260}, {507, 253}, {580, 679, 587, 590, 271, 281, 602}, {272, 585, 592, 356}, {280, 282}, {300, 734, 287}, {304, 377, 360}, {336, 305}, {489, 692, 308, 605}, {313, 338, 357}, {324, 647, 563, 755, 601, 698, 701}, {328, 350, 375}, {737, 733, 616, 332, 621, 688, 625, 756, 568, 730, 635, 509, 511}, {523, 346, 571, 655}, {347, 358}, {362, 541, 349}, {355, 387}, {368, 384}, {745, 386}, {394, 739}, {395, 478}, {412, 579, 404}, {417, 454, 425, 624, 656, 470, 631, 410, 475}, {484, 446, 415}, {433, 527}, {434, 491, 581}, {499, 462}, {472, 674, 467}, {648, 485}, {488, 536, 573, 567}, {524, 764}, {627, 559}, {569, 711}, {570, 653}, {610, 644, 591}, {669, 598, 671}, {600, 646}, {714, 603}, {609, 749, 695}, {637, 687}, {641, 681}, {726, 694}, {738, 708}, {746, 766}, {753, 759}]


Let’s check by visual inspection if the polymers in the second network
are conenected.

.. image:: images/network_rings_stacked.png

These polymers are all conected by polymers with ring stacking!
