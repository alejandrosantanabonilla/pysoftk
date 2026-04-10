import os
import numpy as np
from rdkit import Chem
import deeptime
from deeptime.decomposition import TICA
import hdbscan
import MDAnalysis as mda
from scipy.spatial.distance import cdist
from tqdm.auto import tqdm
from functools import wraps

# Adjusted relative imports for the pol_analysis directory structure
from .tools.utils_mda import MDA_input
from .tools.utils_tools import distance_matrix_calc

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        import time
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        total_time = time.perf_counter() - start_time
        print(f'Function {func.__name__} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper

class AutomatedKineticClustering(MDA_input):
    """
    Automated Kinetic Clustering pipeline for single molecules.
    Uses RDKit for topological backbone identification, Contact Maps for 
    sparse feature extraction, and tICA/HDBSCAN for metastable state clustering.
    """

    def __init__(self, tpr_file, xtc_file):
        super().__init__(tpr_file, xtc_file)
        self.u = self.get_mda_universe()
        
        # Safely extract time step if available, default to 1.0 ps if missing
        self.dt = getattr(self.u.trajectory, 'dt', 1.0) 
        
        self.core_atoms = None
        self.feature_matrix = None
        self.embedding = None
        self.labels = None

    @timeit
    def define_backbone_by_centrality(self, percentile=75):
        """
        Uses RDKit to calculate Closeness Centrality via the topological distance matrix.
        Keeps only the atoms above the given percentile as the 'core backbone'.
        """
        print("Converting MDAnalysis topology to RDKit molecule...")
        try:
            # Native MDA to RDKit conversion preserves atom indices
            rdkit_mol = self.u.atoms.convert_to("RDKIT")
        except Exception as e:
            raise RuntimeError(f"Failed to convert to RDKit. Ensure your topology has bond information. Error: {e}")
        
        print("Calculating RDKit topological distance matrix (C++ optimized)...")
        # Returns a 2D matrix of the shortest path (in number of bonds) between all atom pairs
        dist_matrix = Chem.GetDistanceMatrix(rdkit_mol)
        
        path_sums = np.sum(dist_matrix, axis=1)
        
        # If the user wants the top 25% most central atoms (percentile=75),
        # we need the 25% of atoms with the LOWEST path sum scores.
        threshold_percentile = 100 - percentile
        threshold_value = np.percentile(path_sums, threshold_percentile)
        
        core_indices = np.where(path_sums <= threshold_value)[0]
        self.core_atoms = self.u.atoms[core_indices]
        
        print(f"Reduced {len(self.u.atoms)} total atoms down to a core backbone of {len(self.core_atoms)} atoms (Top {100-percentile}%).")
        return self.core_atoms

    @timeit
    def extract_contact_maps(self, r_cutoff=8.0, step=1):
        """
        Calculates distances for the core atoms and applies a spatial cutoff 
        to create a binary contact matrix for each frame.
        """
        if self.core_atoms is None:
            raise ValueError("Run define_backbone_by_centrality first.")

        print(f"Extracting contact maps with r_cutoff = {r_cutoff} Å...")
        feature_list = []
        
        for ts in tqdm(self.u.trajectory[::step], desc="Building Contact Maps"):
            dist_m = distance_matrix_calc(self.u, self.core_atoms.positions, self.core_atoms.positions)
            
            # Extract upper triangle (exclude diagonal self-distances)
            tri_idx = np.triu_indices(len(self.core_atoms), k=1)
            distances = dist_m[tri_idx]
            
            # Apply cutoff to create a binary matrix (0 or 1)
            contacts = (distances <= r_cutoff).astype(float)
            feature_list.append(contacts)
            
        self.feature_matrix = np.array(feature_list)
        return self.feature_matrix

    @timeit
    def run_tica_clustering(self, lag_time_frames=10, n_components=2, min_cs=15):
        """
        Runs tICA on the binary contact maps and clusters the kinetic space using HDBSCAN.
        """
        if self.feature_matrix is None:
            raise ValueError("Run extract_contact_maps first.")

        physical_lag = lag_time_frames * self.dt
        print(f"Running tICA. Lag time: {lag_time_frames} frames ({physical_lag} ps)...")
        
        tica = TICA(lagtime=lag_time_frames, dim=n_components)
        tica_model = tica.fit(self.feature_matrix).fetch_model()
        self.embedding = tica_model.transform(self.feature_matrix)
        
        print("Clustering kinetic states with HDBSCAN...")
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cs)
        self.labels = clusterer.fit_predict(self.embedding)
        
        n_states = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise = list(self.labels).count(-1)
        print(f"Identified {n_states} distinct metastable states.")
        print(f"Frames in transition states (noise): {n_noise} ({100 * n_noise / len(self.labels):.2f}%)")
        
        return self.embedding, self.labels

    @timeit
    def export_representative_states(self, output_file="metastable_states.pdb", selection="all"):
        """
        Finds the geometric center of each state in the tICA space, locates the 
        closest real frame, and saves them all into a single multi-frame trajectory.
        """
        if self.embedding is None or self.labels is None:
            raise ValueError("Run run_tica_clustering first.")

        unique_labels = set(self.labels)
        sel = self.u.select_atoms(selection)
        
        with mda.Writer(output_file, sel.n_atoms) as W:
            for cluster_id in sorted(unique_labels):
                if cluster_id == -1:
                    continue  # Skip transition states
                
                cluster_indices = np.where(self.labels == cluster_id)[0]
                cluster_points = self.embedding[cluster_indices]
                
                cluster_center = np.mean(cluster_points, axis=0).reshape(1, -1)
                
                distances = cdist(cluster_points, cluster_center)
                closest_local_idx = np.argmin(distances)
                closest_global_idx = cluster_indices[closest_local_idx]
                
                self.u.trajectory[closest_global_idx]
                W.write(sel)
                
                print(f"State {cluster_id} -> Saved Frame {closest_global_idx} to {output_file}")
                
        print("Export complete.")
