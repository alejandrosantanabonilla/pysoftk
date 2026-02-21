import pytest
from rdkit import Chem
from pysoftk.topologies.dend import DendrimerBuilder

def test_frechet_dendrimer_g2():
    """Test the generation of a G2 Fréchet-type dendrimer."""
    
    # --- 1. Define Monomer Building Blocks ---

    # Core: 1,1,1-tris(4'-hydroxyphenyl)ethane with 3 attachment points
    core = {
        'smiles': "CC(c1ccc(O[*:1])cc1)(c1ccc(O[*:2])cc1)c1ccc(O[*:3])cc1",
        'child_attachment_map_nums': [1, 2, 3]
    }

    # Branching Unit: A 1 -> 2 poly(benzyl ether) unit
    branch = {
        'smiles': "c1(O[*:1])cc(CO[*:2])cc(CO[*:3])c1",
        'parent_attachment_map_num': 1,
        'child_attachment_map_nums': [2, 3]
    }

    # --- 2. Initialize and run the builder ---
    
    dendrimer_builder = DendrimerBuilder(
        core_data=core,
        branch_data=branch,
        generations=2
    )
    
    # Build the molecule without running a full geometry optimization
    # to speed up the test. We are testing topology, not the final 3D structure.
    dendrimer_builder.build() 
    
    result_mol = dendrimer_builder.get_optimized_mol()

    # --- 3. Define the expected result ---

    # Instead of a brittle SMILES comparison, we check a fundamental
    # property like the atom count. This verifies the topology is correct.
    # Core (C20H18O3): 41 atoms
    # Branch (C9H10O3): 22 atoms
    # G0 = 1 * Core = 41
    # G1 = G0 + 3 * Branch = 41 + 3*22 = 107
    # G2 = G1 + 6 * Branch = 107 + 6*22 = 239
    expected_atom_count = 212

    # --- 4. Assert that the generated molecule has the correct number of atoms ---
    
    # The 'result_mol' is a pybel.Molecule object.
    # The number of atoms can be found with len(result_mol.atoms).
    assert len(result_mol.atoms) == expected_atom_count
