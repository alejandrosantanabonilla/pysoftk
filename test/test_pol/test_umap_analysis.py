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
from pysoftk.pol_analysis.umap_analysis import umap_analysis

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_umap_anaysis(rootdir):

    cwd=os.path.join(rootdir, 'data/data_umap')

    all_pos=os.path.join(rootdir, 'data/data_umap/all_pos.npy')

    cyclic_x_embed_file=os.path.join(rootdir, 'data/data_umap/cyclic_x_embedding.npy')
    cyclic_x_embed=np.load(cyclic_x_embed_file, allow_pickle=True)

    cyclic_y_pred_file=os.path.join(rootdir, 'data/data_umap/cyclic_y_pred.npy')
    cyclic_y_pred=np.load(cyclic_y_pred_file, allow_pickle=True)

    cyclic_average_rep_file=os.path.join(rootdir, 'data/data_umap/cyclic_av_rep.npy')
    cyclic_average_rep=np.load(cyclic_average_rep_file, allow_pickle=True)

    X_embedded = umap_analysis(all_pos).umap_run(8, cwd)

    y_pred = umap_analysis(all_pos).hdbscan_cluster(X_embedded , 34, 1.1, cwd)

    data_folder = os.path.join(rootdir, 'data/data_umap')
    png_files = [file for file in os.listdir(data_folder) if file.lower().endswith('.png')]
    expected_length = 2
    assert len(png_files) == expected_length










 
