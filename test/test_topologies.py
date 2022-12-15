import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

from pysoftk.topologies.diblock import *
from pysoftk.topologies.ring import *
from pysoftk.topologies.branched import *
from pysoftk.topologies.ranpol import *

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

@pytest.mark.parametrize("mol,expected", testdata)
def test_ring(mol,expected):
    
  mol=Chem.MolFromSmiles(str(mol))
  AllChem.EmbedMolecule(mol)
  
  hom=Rn(mol,'Br').pol_ring(10)
  Fmt(hom).xyz_print("ring.xyz")

  b=Fld().seek_files("xyz")
  assert len(b) == expected
  os.remove("ring.xyz")

def test_branched_1():
    core=Chem.MolFromSmiles('BrN(Br)CCN(Br)Br')
    prob=Chem.MolFromSmiles('[C@H](CCl)(OBr)C')

    # Problematic stereochemistry
    Chem.MolToMolFile(prob,"mol.mol")
    arm=Chem.MolFromMolFile('./mol.mol')

    final = Bd(core, arm, "Br").branched_polymer()
    Fmt(final).xyz_print("bran1.xyz")

    os.remove("mol.mol")
    b=Fld().seek_files("xyz")

    #assert len(b) == 1
    os.remove("bran1.xyz")

def test_branched_2():

    core=Chem.MolFromSmiles('BrN(Br)CCN(Br)Br')
    arm=Chem.MolFromSmiles('BrOCCCl')

    final = Bd(core, arm, "Br").branched_polymer()
    Fmt(final).xyz_print("bran2.xyz")

    b=Fld().seek_files("xyz")

    assert len(b) == 1
    os.remove("bran2.xyz")

def test_ran_pol():
    mol_2=Chem.MolFromSmiles('BrCOCBr')
    mol_4=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
    mol_5=Chem.MolFromSmiles('BrCCBr')

    dia = Rnp(mol_2, mol_4,"Br").random_ab_copolymer(10, 0.4,10)
    Fmt(dia).xyz_print("dia.xyz")
  
    tri = Rnp(mol_2, mol_4,"Br").random_abc_copolymer(mol_2, 15, 0.4, 0.6, 10)
    Fmt(tri).xyz_print("tri.xyz")

    b=Fld().seek_files("xyz")

    assert len(b) == 2
    
    os.remove("dia.xyz")
    os.remove("tri.xyz")
