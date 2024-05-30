import pytest
import os
from functools import wraps
import numpy as np
import filecmp
import pandas as pd

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.ring_ring import RSA
from pysoftk.folder_manager.folder_creator import *

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_RSA(rootdir):

    from pandas.testing import assert_frame_equal

    try:
        topology=os.path.join(rootdir, 'data/trajectory_resids.pdb')
        trajectory=os.path.join(rootdir, 'data/step5_nowater.xtc')

        rsa_og = os.path.join(rootdir, 'data/rsa_small.parquet')
        rsa_og_df= pd.read_parquet(rsa_og)
        results='data/rsa.parquet'
        ang_c=30
        dist_c=5

        rsa=RSA(topology, trajectory).stacking_analysis(dist_c, ang_c, 0, 20, 2, results)
        rsa_df=pd.read_parquet(results)
    
        assert_frame_equal(rsa_df, rsa_og_df)

    except Exception:
        assert False

    finally:
        a=Fld().seek_files("pdb")
        for i in range(len(a)):
            os.unlink(a[i])
  
