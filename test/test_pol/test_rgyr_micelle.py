import pytest
import os
from functools import wraps
import numpy as np 

import pytest
import os
from functools import wraps
import numpy as np

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.rgyr_micelle import rgyr

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))




def test_contacts_cyclic(rootdir):

    cyclic_rgyr_file=os.path.join(rootdir, 'data/cyclic_rgyr.npy')
    cyclic_rgyr=np.load(cyclic_rgyr_file, allow_pickle=True)

    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')

    cyclic_resids = os.path.join(rootdir, 'data/results_cyclic.parquet')

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(cyclic_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)



    MA_oxygen_names = ['O010', 'O012',  'O019', 'O017', 'O00X',  'O00V','O00S', 'O00Z', 'O01I',  'O01G',  'O01E', 'O01C']

    rgyr_micelle_whole = rgyr(topology, trajectory).running_rgyr(resids, atom_pos, 0, 10001, 1)

    plot_rgyr_whole = np.array(list(rgyr_micelle_whole))


    assert np.allclose(plot_rgyr_whole, cyclic_rgyr)



def test_contacts_branched(rootdir):


    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')

    branched_resids = os.path.join(rootdir, 'data/results_branched.parquet')


    branched_rgyr_file=os.path.join(rootdir, 'data/branched_rgyr.npy')
    branched_rgyr_core_file=os.path.join(rootdir, 'data/branched_rgyr_core.npy')
    branched_rgyr=np.load(branched_rgyr_file, allow_pickle=True)
    branched_rgyr_core=np.load(branched_rgyr_core_file, allow_pickle=True)


    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(branched_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)


    MA_names =['O00W', 'O00U', 'O00T', 'O00Y', 'O011', 'O013', 'O01A', 'O018','O01H', 'O01F', 'O01D', 'O01B', 
            'C024', 'C022', 'C020', 'C02A', 'C01X', 'C01T', 'C01P', 'C01L', 'C02Y', 'C02U', 'C02Q', 'C02M', 
          'C023', 'C021', 'C029', 'C02B', 'C01U', 'C01Q', 'C01M', 'C01K', 'C02Z', 'C02V', 'C02R', 'C02N'] 

    rgyr_micelle_whole = rgyr(topology, trajectory).running_rgyr(resids, atom_pos, 0, 10001, 1)
    rgyr_micelle_core = rgyr(topology, trajectory).running_rgyr(resids, atom_pos, 0, 10001, 1, MA_names)

    plot_rgyr_whole = np.array(list(rgyr_micelle_whole))
    plot_rgyr_core = np.array(list(rgyr_micelle_core))

    assert np.allclose(plot_rgyr_whole, branched_rgyr)
    assert np.allclose(plot_rgyr_core, branched_rgyr_core)





 
