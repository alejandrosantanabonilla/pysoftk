import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

from pysoftk.topologies.diblock import *
from pysoftk.topologies.ring import *
from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *

# SMILES topologies polymers testdata
testdata=[
    ('C[Si](OBr)(C)Br',1),
    ('c1cc(ccc1Br)Br',1),
    ('c1(ccc(o1)Br)Br',1),
]


def test_diblock():
    mol_1=Chem.MolFromSmiles('BrCOCBr')
    mol_2=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')
    
    a=Db(mol_1, mol_2, 'Br').diblock_copolymer(5, 2, 10)
    Fmt(a).xyz_print("diblock_monomer.xyz")
    
    b=Fld().seek_files("xyz")
    assert len(b) == 1
    os.remove("diblock_monomer.xyz")

def test_pattern():
    
    mols=[Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'), Chem.MolFromSmiles('BrCCBr')]
    string="ABBAAB"

    a=Pt(string, mols, "Br").pattern_block_poly()
    Fmt(a).xyz_print("patterned_polymer.xyz")

    b=Fld().seek_files("xyz")
    assert len(b) == 1
    os.remove("patterned_polymer.xyz")

@pytest.mark.parametrize("mol,expected", testdata)
def test_ring(mol,expected):
    
  mol=Chem.MolFromSmiles(str(mol))
  AllChem.EmbedMolecule(mol)
  
  hom=Rn(mol,'Br').pol_ring(10)
  Fmt(hom).xyz_print("ring.xyz")

  b=Fld().seek_files("xyz")
  assert len(b) == expected
  os.remove("ring.xyz")
