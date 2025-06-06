import pytest
import os
from functools import wraps
import numpy as np
import filecmp
from openbabel import pybel as pb

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_large_resid_branched(rootdir):
    testdata1 = [([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]),
             ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]),
             ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 18, 19, 20])]

    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')

    branched_resids = os.path.join(rootdir, 'data/results_branched.parquet')

    result = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(branched_resids)

    for i in range(len(testdata1)):
        
        arr = np.array(testdata1[i])  # Convert list to NumPy array
        assert np.array_equal(result[i], arr)


def test_large_resid_cyclic(rootdir):
    testdata2 = [([1, 2, 3, 4, 6, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19]),
       ([2, 4, 6, 7, 13, 14, 15, 17, 18, 19, 20]),
       ([2, 4, 6, 7, 13, 14, 15, 17, 18, 19, 20])]

    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')
    cyclic_resids = os.path.join(rootdir, 'data/results_cyclic.parquet')
    result = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(cyclic_resids)

    for i in range(len(testdata2)):  
        arr = np.array(testdata2[i])  # Convert list to NumPy array
        assert np.array_equal(result[i], arr)


def test_atom_pos_cyclic(rootdir):
    all_pos_cyclic_file=os.path.join(rootdir, 'data/all_atom_pos_cyclic.npy')
    all_pos_cyclic=np.load(all_pos_cyclic_file, allow_pickle=True)
    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')

    cyclic_resids = os.path.join(rootdir, 'data/results_cyclic.parquet')

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(cyclic_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)

    for i in range(len(all_pos_cyclic)):  
        arr = all_pos_cyclic[i] # Convert list to NumPy array
        assert np.array_equal(atom_pos[i][1], arr)


def test_atom_pos_branched(rootdir):
    all_pos_branched_file=os.path.join(rootdir, 'data/all_atom_pos_branched.npy')
    all_pos_branched=np.load(all_pos_branched_file, allow_pickle=True)

    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')

    branched_resids = os.path.join(rootdir, 'data/results_branched.parquet')

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(branched_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)

    for i in range(len(all_pos_branched)):   
        arr = all_pos_branched[i] # Convert list to NumPy array
        assert np.array_equal(atom_pos[i][1], arr)







 
