import pytest
import os
import json
import numpy as np
from umap import UMAP

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.umap_analysis import umap_analysis

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_umap_anaysis(rootdir):

    cwd=os.path.join(rootdir, 'data/data_umap')

    all_pos=os.path.join(rootdir, 'data/data_umap/all_pos.npy')

    # Load UMAP parameters from JSON for X_embedded
    umap_params_file = os.path.join(rootdir, 'data/data_umap/umap_params.json')
    with open(umap_params_file, 'r') as f:
        umap_params = json.load(f)
    
    # Recreate the UMAP object for X_embedded
    umap_object = UMAP(**umap_params) 

    # Assuming you have the data for fitting UMAP stored in 'data_for_umap.npy'
    data_for_umap_file = os.path.join(rootdir, 'data/data_umap/data_for_umap.npy')
    data_for_umap = np.load(data_for_umap_file, allow_pickle=True)

    # Generate X_embedded 
    X_embedded = umap_object.fit_transform(data_for_umap) 


    # Load UMAP parameters from JSON for cyclic_y_pred (same as before)
    with open(umap_params_file, 'r') as f:
        umap_params = json.load(f)
    
    # Recreate the UMAP object for cyclic_y_pred
    cyclic_y_pred = UMAP(**umap_params)

    cyclic_average_rep_file=os.path.join(rootdir, 'data/data_umap/cyclic_av_rep.npy')
    cyclic_average_rep=np.load(cyclic_average_rep_file, allow_pickle=True)

    X_embedded = umap_analysis(all_pos).umap_run(8, cwd)

    y_pred = umap_analysis(all_pos).hdbscan_cluster(X_embedded , 34, 1.1, cwd)

    data_folder = os.path.join(rootdir, 'data/data_umap')
    png_files = [file for file in os.listdir(data_folder) if file.lower().endswith('.png')]
    expected_length = 2
    assert len(png_files) == expected_length

