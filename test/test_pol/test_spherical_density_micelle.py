import pytest
import os
from functools import wraps
import numpy as np

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.spherical_density import spherical_density

@pytest.fixture

def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_sph_density_triblock(rootdir):

    topology=os.path.join(rootdir, 'data/spherical_density.tpr')
    trajectory=os.path.join(rootdir, 'data/spherical_density.xtc')

    resids = os.path.join(rootdir, 'data/spherical.parquet')


    spherical_file=os.path.join(rootdir, 'data/density_PMA.npy')
    spherical_file_f=np.load(spherical_file, allow_pickle=True)


    MA_names_total =['C00A', 'C009', 'C008', 'C007', 'C006', 'C005', 'C004', 'C003', 'C002', 'C001', 'C000', 'C00L', 'C00O', 'C00P', 'C00Q', 'C00R', 'C00S', 'C00T', 
           'C00W', 'C00X', 'C010', 'C011', 'C014', 'C015','O001', 'O002', 'O005', 'O007', 'O008', 'O00A', 'O00O', 'O00L', 'O00C', 'O00F', 'O00G', 'O00J', 
             'O00K', 'C01F', 'C017', 'C013', 'O00H', 'C012', 'C00V', 'O00D', 'C00U', 'C00Z', 
                 'O00E', 'C00Y', 'C019', 'O00M', 'C018',
                 'C01B', 'O00N', 'C01A', 'C00N', 'O00B', 'C00M', 'C00J', 'O009', 'C00K', 'C00H',
                 'O006', 'C00I', 'C00G', 'O004', 'C00F',
                 'C00E', 'O003', 'C00D', 'C00C', 'O000', 'C00B']

    resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(resids)


    atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 0, 1000, 2)


    spherical_density_whole, binned_space = spherical_density(topology, trajectory).run_density_calc('resid ', resids, atom_pos, MA_names_total)

    #spherical_density_whole_f = np.array(list(spherical_density_whole))


    assert np.array_equal(spherical_density_whole, spherical_file_f)
