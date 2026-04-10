import pytest
import os
import pandas as pd
from pandas.testing import assert_frame_equal

from pysoftk.pol_analysis.kinetic_clustering import AutomatedKineticClustering
from pysoftk.folder_manager.folder_creator import *

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_KineticClustering(rootdir):
    try:
        # 1. Define paths to test data (Alanine Dipeptide)
        topology = os.path.join(rootdir, 'data/ala_dip_topology.pdb')
        trajectory = os.path.join(rootdir, 'data/xtb.xyz')
        
        # 2. Path to the Reference "Ground Truth" Data
        ref_labels_file = os.path.join(rootdir, 'data/kinetic_labels_og.parquet')
        ref_labels_df = pd.read_parquet(ref_labels_file)

        # 3. Run the Pipeline
        app = AutomatedKineticClustering(topology, trajectory)
        
        app.define_backbone_by_centrality(percentile=60)
        app.extract_contact_maps(r_cutoff=4.5)
        
        # Run tICA and HDBSCAN
        embedding, labels = app.run_tica_clustering(lag_time_frames=5, n_components=2, min_cs=10)
        
        # Output the states (to test the export function works without crashing)
        output_pdb = os.path.join(rootdir, 'data/test_alanine_states.pdb')
        app.export_representative_states(output_pdb)

        # 4. Convert the output labels to a DataFrame for comparison
        current_labels_df = pd.DataFrame({'cluster_labels': labels})

        # 5. Assert: Does the current run exactly match the reference run?
        assert_frame_equal(current_labels_df, ref_labels_df)

    except Exception as e:
        # Printing the exception helps with debugging if the test fails
        print(f"\nTest failed with exception: {e}")
        assert False

    finally:
        # Cleanup any generated PDB files using your folder manager tool
        try:
            a = Fld().seek_files("pdb")
            for i in range(len(a)):
                # Be careful not to delete the topology file! 
                if "test_alanine_states.pdb" in a[i]:
                    os.unlink(a[i])
        except Exception:
            pass
