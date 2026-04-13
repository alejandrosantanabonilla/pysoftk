import os
import time
from functools import wraps
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances 
from scipy.spatial import cKDTree
from tqdm import tqdm
import itertools
import concurrent.futures

from numba import njit
from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from rdkit import Chem
import MDAnalysis.transformations as trans

import warnings
warnings.filterwarnings('ignore')

# =====================================================================
# NUMBA MATH - nogil=True is the key to Threading speed
# =====================================================================

@njit(nogil=True, fastmath=True, cache=True)
def fast_svd(coordinates):
    """
    Compute the normal vector of a ring plane using Singular Value Decomposition.

    Parameters
    ----------
    coordinates : ndarray of shape (N, 3)
        The Cartesian coordinates of the atoms in the ring.

    Returns
    -------
    ndarray of shape (3,)
        The normal vector (third principal component) of the ring plane.
    """
    mean = np.zeros(3)
    n_atoms = coordinates.shape[0]
    for i in range(n_atoms):
        mean += coordinates[i]
    mean /= n_atoms
    centered = coordinates - mean
    u_, s, vh = np.linalg.svd(centered)
    return np.ascontiguousarray(vh[2, :])

@njit(nogil=True, fastmath=True, cache=True)
def fast_calc_angle(u_norm, u_norm2):
    """
    Calculate the acute angle between two normal vectors.

    Parameters
    ----------
    u_norm : ndarray of shape (3,)
        Normal vector of the first ring.
    u_norm2 : ndarray of shape (3,)
        Normal vector of the second ring.

    Returns
    -------
    float
        The angle between the planes in degrees.
    """
    dot = np.dot(u_norm, u_norm2)
    norm_prod = np.linalg.norm(u_norm) * np.linalg.norm(u_norm2)
    cos_theta = dot / norm_prod
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return np.abs(np.arccos(cos_theta) * (180.0 / np.pi))

# =====================================================================
# THREAD-SAFE WORKER
# =====================================================================

def evaluate_dimer_pair(ring_comb, positions, box_dims, cut_off, ang_cutoff, pol_idx):
    """
    Worker function to check if any pair of rings between two polymers are stacking.

    Parameters
    ----------
    ring_comb : list of tuples
        Combinations of ring atom indices between polymer A and polymer B.
    positions : ndarray of shape (N, 3)
        Global atom positions from the trajectory.
    box_dims : ndarray of shape (6,)
        Simulation box dimensions for PBC handling.
    cut_off : float
        Distance threshold for ring-ring proximity.
    ang_cutoff : float
        Angle threshold (in degrees) for face-to-face stacking.
    pol_idx : tuple
        The Residue IDs of the two polymers being compared.

    Returns
    -------
    tuple or None
        Returns `pol_idx` if a stacking event is found, otherwise None.
    """
    for ringA_idx, ringB_idx in ring_comb:
        posA = positions[ringA_idx]
        posB = positions[ringB_idx]
        
        dist_mat = distances.distance_array(posA, posB, box=box_dims)
        if np.min(dist_mat) < cut_off:
            u_norm = fast_svd(posA)
            u_norm2 = fast_svd(posB)
            angle = fast_calc_angle(u_norm, u_norm2)
            if (angle < ang_cutoff or angle > 180 - ang_cutoff):
                return pol_idx
    return None

# =====================================================================
# MAIN CLASS
# =====================================================================

def timeit(func):
    """Decorator to measure and print function execution time."""
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        print(f'Function {func.__name__} Took {end_time - start_time:.4f} seconds')
        return result
    return timeit_wrapper

class RSA(MDA_input):
    """
    Ring Stacking Analysis (RSA) for polymer systems.

    Inherits from MDA_input to handle MDAnalysis universe initialization.
    """
    def __init__(self, tpr_file, xtc_file):
        super().__init__(tpr_file, xtc_file)
        self.ring_dict = {}

    def _build_global_rings(self, u):
        """
        Identify rings in polymers using RDKit and map them to global MDAnalysis indices.

        Parameters
        ----------
        u : MDAnalysis.Universe
            The universe object containing the molecular system.
        """
        self.ring_dict = {}
        processed_templates = {}
        
        for frag in u.atoms.fragments:
            sig = (tuple(frag.resnames), tuple(frag.names))
            if sig not in processed_templates:
                try:
                    mol = frag.convert_to('RDKIT')
                except:
                    import tempfile
                    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                        frag.write(tmp.name)
                    mol = Chem.MolFromPDBFile(tmp.name, removeHs=False, sanitize=True)
                    os.remove(tmp.name)

                rings = []
                if mol:
                    ring_info = mol.GetRingInfo()
                    for cycle in ring_info.AtomRings():
                        if len(cycle) >= 4:
                            rings.append(cycle)
                processed_templates[sig] = rings

            for local_ring in processed_templates[sig]:
                global_indices = frag.indices[list(local_ring)]
                resid = u.atoms[global_indices[0]].resid
                if resid not in self.ring_dict:
                    self.ring_dict[resid] = []
                self.ring_dict[resid].append(np.sort(global_indices))

    def load_or_build_rings(self, u, cache_file="ring_cache.dill"):
        """
        Load ring data from a cache file or build it if the cache does not exist.

        Parameters
        ----------
        u : MDAnalysis.Universe
            The universe object.
        cache_file : str, optional
            Path to the dill cache file, by default "ring_cache.dill".
        """
        import dill
        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as f:
                self.ring_dict = dill.load(f)
        else:
            self._build_global_rings(u)
            with open(cache_file, 'wb') as f:
                dill.dump(self.ring_dict, f)

    def pol_cutoff(self, u, cutoff):
        """
        Filter polymer pairs that are spatially close using a cKDTree.

        Parameters
        ----------
        u : MDAnalysis.Universe
            The universe object at the current frame.
        cutoff : float
            Distance cutoff for the initial spatial search.

        Returns
        -------
        indices : list
            List of indices (range of valid pairs).
        pairs : list of tuples
            List of polymer Residue ID pairs that are proximity candidates.
        """
        all_ring_indices = [idx for res_rings in self.ring_dict.values() for r in res_rings for idx in r]
        ring_atom_indices = np.unique(all_ring_indices)
        if ring_atom_indices.size == 0: return [], []
        
        ring_atoms = u.atoms[ring_atom_indices]
        box = u.dimensions[:3] if u.dimensions is not None else None
        
        if box is not None:
            wrapped = np.asarray(ring_atoms.positions, dtype=np.float64) % box
            for i in range(3):
                wrapped[:, i] = np.clip(wrapped[:, i], 0.0, box[i] - 1e-5)
            tree = cKDTree(wrapped, boxsize=box)
        else:
            tree = cKDTree(ring_atoms.positions)
            
        pairs = tree.query_pairs(float(cutoff), output_type='ndarray')
        if pairs.size == 0: return [], []
        
        resids = ring_atoms.resids[pairs]
        valid = resids[resids[:,0] != resids[:,1]]
        if valid.size == 0: return [], []
        
        unique_pairs = np.unique(np.sort(valid, axis=1), axis=0)
        return list(range(len(unique_pairs))), [tuple(x) for x in unique_pairs]

    def find_several_rings_stacked(self, rings_df_name):
        """
        Extract connected components (clusters) of polymers from the stacking results.

        Parameters
        ----------
        rings_df_name : str
            Path to the parquet file generated by `stacking_analysis`.

        Returns
        -------
        list of lists
            For each frame, a list of polymer clusters (each cluster is a list of Residue IDs).
        """
        df = pd.read_parquet(rings_df_name)
        out = []
        for i in range(len(df)):
            adj = {}
            raw_pairs = df['pol_resid'].iloc[i]

            for p in raw_pairs:
                u_n, v_n = p[0], p[1]
                adj.setdefault(u_n, set()).add(v_n)
                adj.setdefault(v_n, set()).add(u_n)
            
            seen, comps = set(), []
            for node in adj:
                if node not in seen:
                    c, q = set(), [node]
                    while q:
                        curr = q.pop(0)
                        if curr not in c:
                            c.add(curr); seen.add(curr)
                            q.extend(adj[curr] - c)
                    comps.append(list(c))
            out.append(comps)
        return out

    @timeit
    def stacking_analysis(self, cut_off, ang_cut, start, stop, step, output_name="output.parquet", write_pdb=False):
        """
        Perform high-throughput ring stacking analysis across trajectory frames.

        Parameters
        ----------
        cut_off : float
            Distance threshold for stacking.
        ang_cut : float
            Angle threshold (acute) for face-to-face orientation.
        start : int
            Starting frame index.
        stop : int
            Stopping frame index.
        step : int
            Frame stride.
        output_name : str, optional
            Name of the output parquet file, by default "output.parquet".
        write_pdb : bool, optional
            Whether to write pdb files (legacy flag), by default False.

        Returns
        -------
        pandas.DataFrame
            The results containing polymer stacking pairs per frame.
        """
        u = self.get_mda_universe()
        if not hasattr(u.atoms, 'bonds') or len(u.atoms.bonds) == 0:
            u.atoms.guess_bonds()
        if u.dimensions is not None: 
            u.trajectory.add_transformations(trans.unwrap(u.atoms))
            
        self.load_or_build_rings(u)
        frames = [ts.frame for ts in u.trajectory[int(start):int(stop):int(step)]]
        data = []
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            for f_idx in tqdm(frames, desc="Trajectory Progress"):            
                u.trajectory[f_idx]
                _, dimer_candidates = self.pol_cutoff(u, float(cut_off)) 
                
                pos = u.atoms.positions
                box = u.dimensions
                frame_resids = []
                
                futures = []
                for resA, resB in dimer_candidates:
                    ringsA, ringsB = self.ring_dict.get(resA, []), self.ring_dict.get(resB, [])
                    if ringsA and ringsB:
                        comb = list(itertools.product(ringsA, ringsB))
                        futures.append(executor.submit(evaluate_dimer_pair, comb, pos, box, float(cut_off), float(ang_cut), (resA, resB)))
                
                for future in concurrent.futures.as_completed(futures):
                    pol_idx = future.result()
                    if pol_idx:
                        frame_resids.append(list(pol_idx))
                
                data.append([[], frame_resids])
            
        df = pd.DataFrame(data, columns=['atom_index', 'pol_resid'])
        df.to_parquet(output_name, index=False)
        return df
