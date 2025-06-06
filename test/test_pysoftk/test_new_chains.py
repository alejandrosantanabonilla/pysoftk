from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer.linear_polymer import Lpr

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *
from pysoftk.torsional.torsional import *
from pysoftk.folder_manager.folder_creator import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

import os
import pytest

def test_linear_polymer_benzene():
    """Test simple benzene polymerization."""
    smiles = "c1ccccc1{R}"
    replacements = {"R": "c2ccc(cc2){R}"}
    user_final_replacement = ""
    lpr_instance = Lpr(smiles, replacements, max_repetitions=9,final_replacement=user_final_replacement)
    mol = lpr_instance.generate_recursive_smiles()
    res = Chem.MolFromSmiles(mol.write("smi"))

    # Generate the expected SMILES based on the benzene polymerization logic.
    expected_smiles = "c1ccccc1c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccccc1"


    assert Chem.MolToSmiles(res) == Chem.MolToSmiles(Chem.MolFromSmiles(expected_smiles))


def test_linear_polymer_styrene():
    """Test styrene polymerization."""
    smiles2 = "CC(c1ccccc1){R}"
    replacements2 = {"R": "CC(c1ccccc1){R}"}
    final_replacement2 = ""
    lpr_instance2 = Lpr(smiles2, replacements2, max_repetitions=10, final_replacement=final_replacement2) # reduce max_repetitions for faster test.
    result2 = lpr_instance2.generate_recursive_smiles(force_field="MMFF", relax_iterations=50, rot_steps=1) # Reduce relax_iterations for faster test.
    res2 = Chem.MolFromSmiles(result2.write("smi"))

   # Generate the expected SMILES based on the styrene polymerization logic.
    expected_smiles = "CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CCc1ccccc1"

    assert Chem.MolToSmiles(res2) == Chem.MolToSmiles(Chem.MolFromSmiles(expected_smiles))

def test_linear_polymer_ethylene():
    """Test ethylene polymerization."""
    smiles3 = "C=C{R}"
    replacements3 = {"R": "C=C{R}"}
    final_replacement3 = ""
    lpr_instance3 = Lpr(smiles3, replacements3, max_repetitions=20, final_replacement=final_replacement3)
    result3 = lpr_instance3.generate_recursive_smiles(force_field="MMFF", relax_iterations=50, rot_steps=1) # reduce relax_iterations for faster test.
    res3 = Chem.MolFromSmiles(result3.write("smi"))
    compare3 = "C=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=C"

    assert Chem.MolToSmiles(res3) == Chem.MolToSmiles(Chem.MolFromSmiles(compare3))
