import pytest
import os
from functools import wraps
import numpy as np

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.ecc_micelle import ecc

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_ecc_cyclic(rootdir):

    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')

    cyclic_resids = os.path.join(rootdir, 'data/results_cyclic.parquet')


    cyclic_ecc_file=os.path.join(rootdir, 'data/cyclic_ecc.npy')
    cyclic_ecc=np.load(cyclic_ecc_file, allow_pickle=True)

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(cyclic_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)


    ecc_micelle_whole = ecc(topology, trajectory).running_ecc(resids, atom_pos, 0, 10001, 1)

    plot_ecc_whole = np.array(list(ecc_micelle_whole))


    assert np.array_equal(plot_ecc_whole, cyclic_ecc)



def test_ecc_branched(rootdir):
    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')

    branched_resids=os.path.join(rootdir, 'data/results_branched.parquet')
    branched_ecc_file=os.path.join(rootdir, 'data/branched_ecc.npy')

    branched_ecc=np.load(branched_ecc_file, allow_pickle=True)

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(branched_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)


    MA_names =['O00W', 'O00U', 'O00T', 'O00Y', 'O011', 'O013', 'O01A', 'O018','O01H', 'O01F', 'O01D', 'O01B', 
            'C024', 'C022', 'C020', 'C02A', 'C01X', 'C01T', 'C01P', 'C01L', 'C02Y', 'C02U', 'C02Q', 'C02M', 
          'C023', 'C021', 'C029', 'C02B', 'C01U', 'C01Q', 'C01M', 'C01K', 'C02Z', 'C02V', 'C02R', 'C02N'] 

    ecc_micelle_whole = ecc(topology, trajectory).running_ecc(resids, atom_pos, 0, 10001, 1)
    plot_ecc_whole = np.array(list(ecc_micelle_whole))

    assert np.array_equal(plot_ecc_whole, branched_ecc)





 
