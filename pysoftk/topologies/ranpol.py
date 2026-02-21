"""
Optimized Polymer Builder Module.

This module provides the `Rnp` class for the computational synthesis of 
ultra-large random copolymers. It utilizes a highly efficient binary tree 
assembly algorithm for topology construction and integrates OpenBabel 
for realistic 3D spatial geometry optimization.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

# 1. SILENCE WARNINGS
RDLogger.DisableLog('rdApp.*')

# Robust Import for Super Monomer
try:
    from pysoftk.linear_polymer.super_monomer import Sm
except ImportError:
    try:
        from super_monomer import Sm
    except ImportError:
        from .super_monomer import Sm

# Import external OpenBabel utilities
try:
    from pysoftk.tools.utils_ob import *
except ImportError:
    print("Warning: Could not import pysoftk.tools.utils_ob")


class Rnp:
    """
    Optimized Polymer Builder for Ultra-Large Molecules.

    This class constructs A-B and A-B-C random copolymers. It prevents the 
    typical O(N^2) scaling bottlenecks of linear polymer growth by using a 
    binary tree assembly approach. It also delegates 3D coordinate generation 
    to OpenBabel for physically realistic spatial accommodation.

    Attributes
    ----------
    ma : rdkit.Chem.rdchem.Mol
        RDKit Mol object representing monomer A.
    mb : rdkit.Chem.rdchem.Mol
        RDKit Mol object representing monomer B.
    atom : str
        The placeholder atomic symbol (e.g., 'Br', 'At') used as the 
        connection point between monomers.
    """

    __slots__ = ['ma', 'mb', 'atom']

    def __init__(self, ma: Chem.Mol, mb: Chem.Mol, atom: str):
        """
        Initializes the Rnp builder with base monomers and a placeholder atom.

        Parameters
        ----------
        ma : rdkit.Chem.rdchem.Mol
            The first monomer template.
        mb : rdkit.Chem.rdchem.Mol
            The second monomer template.
        atom : str
            The atomic symbol used as the reactive placeholder.
        """
        self.ma = ma
        self.mb = mb
        self.atom = atom

    def _build_sequence_tree(self, mol_list: list) -> Chem.Mol:
        """
        Assembles a list of molecules into a single chain using binary tree reduction.
        

        Instead of adding monomers one by one linearly, this pairs them up and 
        merges them recursively, reducing assembly depth to O(log N).

        Parameters
        ----------
        mol_list : list of rdkit.Chem.rdchem.Mol
            An ordered sequence of monomer RDKit Mol objects.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            A single RDKit Mol object containing the 2D topology of the entire chain.
        """
        current_layer = mol_list
        while len(current_layer) > 1:
            next_layer = []
            
            for i in range(0, len(current_layer) - 1, 2):
                mol_left = current_layer[i]
                mol_right = current_layer[i+1]
                
                # Instant topology construction
                combined = Sm(mol_left, mol_right, str(self.atom)).build_topology_only()
                next_layer.append(combined)

            # Carry over the odd molecule if the layer has an uneven count
            if len(current_layer) % 2 == 1:
                next_layer.append(current_layer[-1])
            current_layer = next_layer
            
        return current_layer[0]

    def _swap_h_topology_only(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Replaces remaining terminal placeholder atoms with Hydrogen purely at the graph level.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The constructed polymer chain containing terminal placeholder atoms.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The modified polymer with placeholder atoms replaced by Hydrogen.
        """
        if not mol:
            return None
            
        pattern = Chem.MolFromSmiles(f"[{self.atom}]")
        replacement = Chem.MolFromSmiles("[H]")
        
        try:
            res = AllChem.ReplaceSubstructs(mol, pattern, replacement, replaceAll=True)
            new_mol = res[0] if res else mol
            Chem.SanitizeMol(new_mol)
            return new_mol
        except Exception:
            return mol

    def _finalize_structure(self, mol: Chem.Mol, relax_iterations: int, 
                            force_field: str, swap_H: bool, 
                            rot_steps: int, ff_thr: float) -> Chem.Mol:
        """
        Handles final topological cleanup and bridges to OpenBabel for 3D generation.
        

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The 2D polymer graph.
        relax_iterations : int
            Number of iterations for the OpenBabel geometry optimization.
        force_field : str
            The force field to apply (e.g., 'MMFF94', 'UFF').
        swap_H : bool
            Whether to replace terminal placeholders with Hydrogen atoms.
        rot_steps : int
            Number of iterations for the weighted rotor search.
        ff_thr : float
            Convergence threshold for the optimization.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The final 3D-optimized polymer as a standard RDKit Mol object.
        """
        final_mol = mol

        # Ensure Hydrogens are present for accurate 3D physics
        if final_mol:
            final_mol = Chem.AddHs(final_mol)

        # Swap Placeholders for H
        if swap_H:
            final_mol = self._swap_h_topology_only(final_mol)
            if final_mol:
                final_mol = Chem.AddHs(final_mol)

        # --- 3D Generation via OpenBabel Bridge ---
        if final_mol:
            try:
                from openbabel import pybel as pb
                
                mol_block_2d = Chem.MolToMolBlock(final_mol)
                ob_mol = pb.readstring("mol", mol_block_2d)
                
                # Map force field names correctly
                ff_lower = force_field.lower() if force_field else "mmff94"
                if ff_lower == "mmff":
                    ff_lower = "mmff94"
                    
                ff_standard = "MMFF94" if ff_lower == "mmff94" else force_field.upper()
                
                print(f"\n[Physics] Generating spatial 3D coordinates via OpenBabel ({ff_lower})...")
                # Step 1: Fast initial spatial placement out of the 2D plane
                ob_mol.make3D(forcefield=ff_lower, steps=50)
                
                # Step 2: Global Optimization (Relaxation + Rotor Search)
                if relax_iterations > 0:
                    print(f"[Physics] Calling external global_opt (iters: {relax_iterations}, rotor_steps: {rot_steps}, thr: {ff_thr})...")
                    ob_mol = global_opt(ob_mol, FF=ff_standard, relax_iterations=relax_iterations, 
                                        rot_steps=rot_steps, ff_thr=ff_thr)
                
                # Bridge back to RDKit
                mol_block_3d = ob_mol.write("mol")
                final_mol_3d = Chem.MolFromMolBlock(mol_block_3d, removeHs=False)
                
                return final_mol_3d if final_mol_3d else final_mol 

            except Exception as e:
                print(f"\n[Warning] OpenBabel Bridge failed ({e}). Falling back to RDKit's random 3D coordinates.")
                try:
                    res = AllChem.EmbedMolecule(final_mol, useRandomCoords=True)
                    if res != -1 and relax_iterations > 0:
                        try:
                            if force_field in ["MMFF", "MMFF94", "mmff94"]:
                                AllChem.MMFFOptimizeMolecule(final_mol, maxIters=relax_iterations)
                            else:
                                AllChem.UFFOptimizeMolecule(final_mol, maxIters=relax_iterations)
                        except Exception:
                            pass 
                except Exception:
                    pass 

        return final_mol

    def random_ab_copolymer(self, len_polymer: int, pA: float, 
                            relax_iterations: int = 100, force_field: str = "MMFF94", 
                            swap_H: bool = True, rot_steps: int = 125, 
                            ff_thr: float = 1.0e-6) -> Chem.Mol:
        """
        Constructs a randomly sequenced two-component (A-B) copolymer.

        Parameters
        ----------
        len_polymer : int
            The total number of monomer units in the desired polymer chain.
        pA : float
            The probability (0.0 to 1.0) of selecting monomer A at any position.
        relax_iterations : int, optional
            Number of iterations for geometry relaxation. Defaults to 100.
        force_field : str, optional
            The force field for 3D optimization. Defaults to 'MMFF94'.
        swap_H : bool, optional
            If True, caps terminal placeholders with Hydrogen. Defaults to True.
        rot_steps : int, optional
            Thoroughness of the OpenBabel rotor search. Defaults to 125.
        ff_thr : float, optional
            Convergence threshold for optimization. Defaults to 1.0e-6.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The synthesized and 3D-optimized A-B copolymer.
        """
        monomers = [self.ma, self.mb]
        probs = [pA, 1.0 - pA]
        
        sequence_indices = np.random.choice([0, 1], size=len_polymer, p=probs)
        mol_sequence = [monomers[i] for i in sequence_indices]

        raw_poly = self._build_sequence_tree(mol_sequence)
        return self._finalize_structure(raw_poly, relax_iterations, force_field, swap_H, rot_steps, ff_thr)

    def random_abc_copolymer(self, mc: Chem.Mol, len_polymer: int, pA: float, pB: float, 
                             relax_iterations: int = 100, force_field: str = "MMFF94", 
                             swap_H: bool = True, rot_steps: int = 125, 
                             ff_thr: float = 1.0e-6) -> Chem.Mol:
        """
        Constructs a randomly sequenced three-component (A-B-C) terpolymer.

        Parameters
        ----------
        mc : rdkit.Chem.rdchem.Mol
            RDKit Mol object representing the third monomer, C.
        len_polymer : int
            The total number of monomer units in the desired polymer chain.
        pA : float
            The probability (0.0 to 1.0) of selecting monomer A.
        pB : float
            The probability (0.0 to 1.0) of selecting monomer B. 
            (Note: Probability of C is calculated as 1.0 - pA - pB).
        relax_iterations : int, optional
            Number of iterations for geometry relaxation. Defaults to 100.
        force_field : str, optional
            The force field for 3D optimization. Defaults to 'MMFF94'.
        swap_H : bool, optional
            If True, caps terminal placeholders with Hydrogen. Defaults to True.
        rot_steps : int, optional
            Thoroughness of the OpenBabel rotor search. Defaults to 125.
        ff_thr : float, optional
            Convergence threshold for optimization. Defaults to 1.0e-6.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The synthesized and 3D-optimized A-B-C terpolymer.
        """
        pC = max(0.0, 1.0 - (pA + pB))
        probs = np.array([pA, pB, pC])
        probs /= probs.sum()

        choices = [self.ma, self.mb, mc]
        indices = np.random.choice([0, 1, 2], size=len_polymer, p=probs)
        mol_sequence = [choices[i] for i in indices]

        raw_poly = self._build_sequence_tree(mol_sequence)
        return self._finalize_structure(raw_poly, relax_iterations, force_field, swap_H, rot_steps, ff_thr)
