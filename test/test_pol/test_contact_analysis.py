import pytest
import os
import numpy as np
from functools import wraps

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.contact_analysis import contacts

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_contacts_cyclic(rootdir):

    contacts_cyclic_file=os.path.join(rootdir, 'data/cyclic_ma_ma_contacts.npy')
    contacts_cyclic=np.load(contacts_cyclic_file, allow_pickle=True)

    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')

    cyclic_resids = os.path.join(rootdir, 'data/results_cyclic.parquet')

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(cyclic_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)

    MA_names = ['C027',  'C023',  'C021',  'C02H',  'C02L',  'C02P',  'C02T', 'C02X',
                'C00U',  'C00R',  'C00P', 'C00L']

    contacts_matrix = contacts(topology, trajectory).run_contacts_calc(resids, atom_pos,
                                                                       MA_names,MA_names, 10)

    contacts_f = np.array(list(contacts_matrix))

    result = []

    for i in range(len(contacts_f)):
        result.append(list(contacts_f[i]))
    
    results = np.array(result)

    contacts_total = np.sum(result, axis=0)

    for i in range(len(contacts_total)):
        arr = np.array(contacts_cyclic[i])  # Convert list to NumPy array
        assert np.array_equal(contacts_total[i], arr)


def test_contacts_branched(rootdir):

    contacts_branched_file=os.path.join(rootdir, 'data/branched_ma_ma_contacts.npy')

    contacts_branched=np.load(contacts_branched_file, allow_pickle=True)


    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')

    branched_resids = os.path.join(rootdir, 'data/results_branched.parquet')

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(branched_resids)
    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 10001, 1)

    MA_names=['C024', 'C022', 'C020', 'C02A', 'C01X', 'C01T', 'C01P', 'C01L', 'C02Y', 'C02U', 'C02Q', 'C02M']

    contacts_matrix = contacts(topology, trajectory).run_contacts_calc(resids, atom_pos, MA_names,MA_names, 10)

    contacts_f = np.array(list(contacts_matrix))

    result = []

    for i in range(len(contacts_f)):
        result.append(list(contacts_f[i]))
    
    results = np.array(result)

    contacts_total = np.sum(result, axis=0)

    for i in range(len(contacts_total)):
        
        arr = np.array(contacts_branched[i])  # Convert list to NumPy array
        assert np.array_equal(contacts_total[i], arr)
