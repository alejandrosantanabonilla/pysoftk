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
testdata=[('C[Si](OBr)(C)Br',1),
    ('c1cc(ccc1Br)Br',1),
    ('c1(ccc(o1)Br)Br',1)]

testdata2=[('c1cc(ccc1Br)Br','c1cc2ccc1-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc-2cc1'),
           ('c1(ccc(o1)Br)Br','c1cc2oc1-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc-2o1')]
    
@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_diblock_1(rootdir):
    est_file = os.path.join(rootdir, 'tmp_folder11')
    os.mkdir(est_file)
    os.chdir(est_file)

    mol_1=Chem.MolFromSmiles('BrCOCBr')
    mol_2=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')
    
    a=Db(mol_1, mol_2, 'Br').diblock_copolymer(5, 2, "MMFF94", 10)
    Fmt(a).xyz_print("diblock_monomer.xyz")
    
    b=Fld().seek_files("xyz")
    assert len(b) == 1

    shutil.rmtree(est_file)

def test_diblock_2(rootdir):
    est_file = os.path.join(rootdir, 'tmp_folder29')
    os.mkdir(est_file)
    os.chdir(est_file)

    mol_1=Chem.MolFromSmiles('BrCOCBr')
    mol_2=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')
    
    a=Db(mol_1, mol_2, 'Br').diblock_copolymer(6, 3, "UFF", 10)
    Fmt(a).xyz_print("diblock_monomer.xyz")
    
    b=Fld().seek_files("xyz")
    assert len(b) == 1

    shutil.rmtree(est_file)
    
@pytest.mark.parametrize("mol,expected", testdata2)
def test_ring(mol, expected, rootdir):
          
  mol=Chem.MolFromSmiles(str(mol))
  AllChem.EmbedMolecule(mol)
  hom=Rn(mol,'Br').pol_ring(10,"MMFF",250)
  new=Chem.MolFromSmiles(str(hom))
  
  assert Chem.MolToSmiles(new) == expected


def test_branched_1(rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder12')
    os.mkdir(est_file)
    os.chdir(est_file)

    core=Chem.MolFromSmiles('BrN(Br)CCN(Br)Br')
    prob=Chem.MolFromSmiles('[C@H](CCl)(OBr)C')

    # Problematic stereochemistry
    Chem.MolToMolFile(prob,"mol.mol")
    arm=Chem.MolFromMolFile('./mol.mol')

    final = Bd(core, arm, "Br").branched_polymer()
    Fmt(final).xyz_print("bran1.xyz")

    b=Fld().seek_files("xyz")

    assert len(b) == 1
    
    shutil.rmtree(est_file)

def test_branched_2(rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder13')
    os.mkdir(est_file)
    os.chdir(est_file)

    core=Chem.MolFromSmiles('BrN(Br)CCN(Br)Br')
    arm=Chem.MolFromSmiles('BrOCCCl')

    final = Bd(core, arm, "Br").branched_polymer()
    Fmt(final).xyz_print("bran2.xyz")

    b=Fld().seek_files("xyz")

    assert len(b) == 1

    shutil.rmtree(est_file)

def test_ran_pol(rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder13')
    os.mkdir(est_file)
    os.chdir(est_file)
    
    mol_2=Chem.MolFromSmiles('BrCOCBr')
    mol_4=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
    mol_5=Chem.MolFromSmiles('BrCCBr')

    dia = Rnp(mol_2, mol_4,"Br").random_ab_copolymer(10, 0.4,force_field="MMFF94",iter_ff=10)
    Fmt(dia).xyz_print("dia.xyz")
  
    tri = Rnp(mol_2, mol_4,"Br").random_abc_copolymer(mol_2, 15, 0.4, 0.6, force_field="MMFF94",iter_ff=10)
    Fmt(tri).xyz_print("tri.xyz")

    b=Fld().seek_files("xyz")

    assert len(b) == 2
    
    shutil.rmtree(est_file)
