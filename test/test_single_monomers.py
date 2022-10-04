import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *

# SMILES asymmetrical polymers testdata
testdata=[
    ('c1cc(sc1CBr)Br',1),
    ('BrOCCBr',1),
    ('[C@@H](CBr)(Br)C(=O)N',1),
    ('C([C@@H](CN)Br)Br',1),
    ('[C@@H](C(=O)OBr)(C)Br',1),
]

# SMILES asymmetrical polymers testdata_2
testdata_2=[
    ('[C@@H](C(=O)OBr)(C)Br',1),
    ('[C@@H](CBr)(Br)C(=O)OC',1),
    ('[C@](CBr)(Br)(C(=O)OC)C',1),
    ('c1c([C@@H](CBr)Br)cccc1',1),
]

@pytest.mark.parametrize("mol,expected", testdata)
def test_single_monomer(mol,expected):
    
    # Original Embedding
    a=Chem.MolFromSmiles(str(mol))
    AllChem.EmbedMolecule(a)

    new=Lp(a,"Br",7,shift=1.25).linear_polymer("MMFF",350)
    Fmt(new).xyz_print("mol_1.xyz")

    b=Fld().seek_files("xyz")
    assert len(b) == expected
    os.remove("mol_1.xyz")

@pytest.mark.parametrize("mol,expected", testdata_2)
def test_single_monomer_2(mol, expected):
    
    a=Chem.MolFromSmiles(str(mol))

    #New_embedding
    ps = rdDistGeom.ETKDGv3()
    AllChem.EmbedMolecule(a,ps)

    new=Lp(a,"Br",7,shift=1.25).linear_polymer("MMFF",350)
    Fmt(new).xyz_print("mol_1.xyz")

    b=Fld().seek_files("xyz")
    
    assert len(b) == expected
    os.remove("mol_1.xyz")

def test_problematic_1():
    # Procedure for molecules with a very complicated stereochemistry
    # One needs to read in SMILES, print into MOL and reload in SMILES
    # MOL-format enables the creation of bonds and cures this issue.
    
    a=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')
    b=Chem.MolToSmiles(a, isomericSmiles=False)
    c=Chem.MolFromSmiles(b)

    ps = rdDistGeom.ETKDGv3()
    AllChem.EmbedMolecule(c,ps)

    new=Lp(c,"Br",7,shift=1.25).linear_polymer("MMFF",350)
    Fmt(new).xyz_print("mol_1.xyz")

    b=Fld().seek_files("xyz")
    
    assert len(b) == 1
    os.remove("mol_1.xyz")
    
def test_problematic_2():
    # Last example of a molecule withcomplicated definition.
    
    a=Chem.MolFromSmiles("C(/C=C(\C)/CBr)Br")

    ps = rdDistGeom.ETKDGv3()
    AllChem.EmbedMolecule(a,ps)

    new=Lp(a,"Br",7,shift=1.25).linear_polymer("MMFF",350)
    Fmt(new).xyz_print("mol_1.xyz")

    b=Fld().seek_files("xyz")
    
    assert len(b) == 1
    os.remove("mol_1.xyz")
    
