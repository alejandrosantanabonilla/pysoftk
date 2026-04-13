import pytest
import os
import pandas as pd
import numpy as np
from pysoftk.pol_analysis.ring_ring import RSA
from pysoftk.folder_manager.folder_creator import *

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_RSA(rootdir):
    # 1. Setup paths
    topology = os.path.join(rootdir, 'data/trajectory_resids.pdb')
    trajectory = os.path.join(rootdir, 'data/step5_nowater.xtc')
    rsa_og_path = os.path.join(rootdir, 'data/rsa_small.parquet')
    
    results = 'data/rsa.parquet'
    cache_file = 'ring_cache.dill'
    
    ang_c = 30
    dist_c = 5

    try:
        # 2. Run the actual analysis
        rsa_calc = RSA(topology, trajectory)
        rsa_calc.stacking_analysis(dist_c, ang_c, 0, 20, 2, results)
        
        # 3. Load results for comparison
        rsa_df = pd.read_parquet(results)
        rsa_og_df = pd.read_parquet(rsa_og_path)

        # 4. Helper to normalize data structure and order
        def normalize_pairs(row_data):
            flat_pairs = []
            def find_pairs(item):
                if isinstance(item, (list, np.ndarray, tuple)):
                    if len(item) == 2 and not isinstance(item[0], (list, np.ndarray, tuple)):
                        flat_pairs.append(tuple(sorted(item)))
                    else:
                        for sub in item:
                            find_pairs(sub)
            find_pairs(row_data)
            return set(flat_pairs)

        # 5. Validation loop
        for i in range(len(rsa_df)):
            new_pairs = normalize_pairs(rsa_df['pol_resid'].iloc[i])
            old_pairs = normalize_pairs(rsa_og_df['pol_resid'].iloc[i])
            
            assert new_pairs == old_pairs, f"Mismatch in frame {i}."

        print("\n✅ Test Passed: Results are consistent.")

    finally:
        # --- THE CLEANUP CREW ---
        # This code runs even if the test fails above
        
        # 1. Erase the ring cache
        if os.path.exists(cache_file):
            os.remove(cache_file)
            print(f"\nSuccessfully erased {cache_file}")

        # 2. Erase the generated parquet
        if os.path.exists(results):
            os.remove(results)

        # 3. Erase any temporary PDB snapshots
        # (Using your folder_creator utility)
        pdbs = Fld().seek_files("pdb")
        for f in pdbs:
            os.unlink(f)
  
