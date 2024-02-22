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
from pysoftk.pol_analysis.contact_analysis import contacts
from pysoftk.pol_analysis.solvation import solvation

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_hydration_cyclic(rootdir):

    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')

    cyclic_resids = os.path.join(rootdir, 'data/results_cyclic.parquet')


    hydration_cyclic_file=os.path.join(rootdir, 'data/cyclic_hydration.npy')
    hydration_cyclic=np.load(hydration_cyclic_file, allow_pickle=True)

    water_topology=os.path.join(rootdir, 'data/cyclic_water_short_trajectory.tpr')
    water_trajectory=os.path.join(rootdir, 'data/cyclic_water_short_trajectory.xtc')


    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(cyclic_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)

    MA_oxygen_names = ['O00A', 'O00D', 'O00C', 'O01G', 'O01F', 'O01D', 'O01A', 'O019', 'O016', 'O015', 'O010', 'O012']

    coord_number = solvation(water_topology, water_trajectory).solvation_calc_run(0, 10001, 1, resids, atom_pos, ['OW'], MA_oxygen_names, 4.5 )
    

    coord_f = np.array(list(coord_number))

    coord_f_r=np.sum(coord_f, axis=0)



    assert set(coord_f_r) == set(hydration_cyclic)


def test_hydration_branched(rootdir):


    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')

    branched_resids = os.path.join(rootdir, 'data/results_branched.parquet')


    hydration_branched_file=os.path.join(rootdir, 'data/branched_hydration.npy')
    hydration_branched=np.load(hydration_branched_file, allow_pickle=True)

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(branched_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)



    MA_oxygen_names = ['O010', 'O012',  'O019', 'O017', 'O00X',  'O00V','O00S', 'O00Z', 'O01I',  'O01G',  'O01E','O01C']

    coord_number = solvation(topology, trajectory).solvation_calc_run(0, 10001, 1, resids, atom_pos, ['OW'], MA_oxygen_names, 4.5 )
    

    coord_f = np.array(list(coord_number))

    coord_f_r=np.sum(coord_f, axis=0)

    assert set(coord_f_r) == set(hydration_branched)
