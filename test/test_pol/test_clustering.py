import pytest
from functools import wraps
import pandas as pd
from pandas.testing import assert_frame_equal
import os

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_clustering_branched(rootdir):

    test_branched = os.path.join(rootdir, 'data/results_branched.parquet')
    testdata1=pd.read_parquet(test_branched)

    atom_names = ['C02B', 'C01K', 'C02N']
    cluster_cutoff = 12

    topology=os.path.join(rootdir, 'data/short_movie_branched.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_branched.xtc')
    results_name=os.path.join(rootdir, 'data/results_branched_2')

    c = SCP(topology, trajectory).spatial_clustering_run(0, 10001, 1, atom_names, cluster_cutoff, results_name)


    df_results = os.path.join(rootdir, 'data/results_branched_2.parquet')
    df = pd.read_parquet(df_results)
    
    # Select the rows for comparison
    expected_values = testdata1.loc[testdata1['time'].isin(df['time'])]

    assert_frame_equal(df, expected_values)

    os.unlink(df_results)



def test_clustering_cyclic(rootdir):

    test_cyclic = os.path.join(rootdir, 'data/results_cyclic.parquet')
    testdata2=pd.read_parquet(test_cyclic)

    atom_names=['C02T']
    cluster_cutoff = 12


    topology=os.path.join(rootdir, 'data/short_movie_cyclic.tpr')
    trajectory=os.path.join(rootdir, 'data/short_movie_cyclic.xtc')
    results_name=os.path.join(rootdir, 'data/results_cyclic_2')

    c = SCP(topology, trajectory).spatial_clustering_run(0, 10001, 1, atom_names, cluster_cutoff, results_name)

    df_results = os.path.join(rootdir, 'data/results_cyclic_2.parquet')
    df = pd.read_parquet(df_results)
    
    # Select the rows for comparison
    expected_values = testdata2.loc[testdata2['time'].isin(df['time'])]

    assert_frame_equal(df, expected_values)#
    
    os.unlink(df_results)
