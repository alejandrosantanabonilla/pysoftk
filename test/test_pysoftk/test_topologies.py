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

from openbabel import openbabel as ob
from openbabel import pybel as pb

# SMILES topologies polymers testdata
testdata=[('C[Si](OBr)(C)Br',1),
    ('c1cc(ccc1Br)Br',1),
    ('c1(ccc(o1)Br)Br',1)]

testdata2=[('c1cc(ccc1Br)Br','c1cc2ccc1-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc-2cc1'),
           ('c1(ccc(o1)Br)Br','c1cc2oc1-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc(o1)-c1ccc-2o1')]

testdata3=[('[C@H](CBr)(OBr)C',
            'C[C@@H](CO[C@@H](C)CCOCCOCCOCCOCCOC)O')]


testdata4=[('[C@H](CCl)(OBr)C',
            '[H]C([H])(N(O[C@@]([H])(C([H])([H])[H])C([H])([H])Cl)O[C@@]([H])(C([H])([H])[H])C([H])([H])Cl)C([H])([H])N(O[C@@]([H])(C([H])([H])[H])C([H])([H])Cl)O[C@@]([H])(C([H])([H])[H])C([H])([H])Cl')]

testdata5=[('BrN(Br)CCN(Br)Br',
            '[H]C([H])(Cl)C([H])([H])ON(OC([H])([H])C([H])([H])Cl)C([H])([H])C([H])([H])N(OC([H])([H])C([H])([H])Cl)OC([H])([H])C([H])([H])Cl')]

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

@pytest.mark.parametrize("comb,expected", testdata3)
def test_diblock_1(comb, expected):
    mol_1=Chem.MolFromSmiles('BrCOCBr')
    mol_2=Chem.MolFromSmiles(comb)
    
    a=Db(mol_1, mol_2, 'Br').diblock_copolymer(5, 2)
    res=Chem.MolFromSmiles(a.write("smi"))

    assert Chem.MolToSmiles(res) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))

def test_diblock_2(rootdir):
    """
    Testing relax_iterations new flag
    """
    est_file = os.path.join(rootdir, 'tmp_folder29')
    os.mkdir(est_file)
    os.chdir(est_file)

    mol_1=Chem.MolFromSmiles('BrCOCBr')
    mol_2=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')
    
    a=Db(mol_1, mol_2, 'Br').diblock_copolymer(6, 3, "UFF", relax_iterations=10)
    Fmt(a).xyz_print("diblock_monomer.xyz")
    
    b=Fld().seek_files("xyz")
    assert len(b) == 1

    shutil.rmtree(est_file)
    
@pytest.mark.parametrize("mol,expected", testdata2)
def test_ring(mol, expected, rootdir):
          
  mol=Chem.MolFromSmiles(str(mol))
  AllChem.EmbedMolecule(mol)
  hom=Rn(mol,'Br').pol_ring(10,"MMFF",relax_iterations=250)
  new=Chem.MolFromSmiles(str(hom))
  
  assert Chem.MolToSmiles(new) == expected

    
@pytest.mark.parametrize("mol,expected", testdata4)
def test_branched_1(mol, expected, rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder12')
    os.mkdir(est_file)
    os.chdir(est_file)

    core=Chem.MolFromSmiles('BrN(Br)CCN(Br)Br')
    prob=Chem.MolFromSmiles(str(mol))

    # Problematic stereochemistry
    Chem.MolToMolFile(prob,"mol.mol")
    arm=Chem.MolFromMolFile('./mol.mol')

    final = Bd(core, arm, "Br").branched_polymer()
    Fmt(final).xyz_print("bran1.xyz")

    b=Fld().seek_files("xyz")
    assert len(b) == 1

    assert Chem.MolToSmiles(final) == expected
    
    shutil.rmtree(est_file)

@pytest.mark.parametrize("mol,expected", testdata5)
def test_branched_2(mol, expected, rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder13')
    os.mkdir(est_file)
    os.chdir(est_file)

    core=Chem.MolFromSmiles(str(mol))
    arm=Chem.MolFromSmiles('BrOCCCl')

    final = Bd(core, arm, "Br").branched_polymer()
    Fmt(final).xyz_print("bran2.xyz")

    b=Fld().seek_files("xyz")
           
    assert len(b) == 1
    assert Chem.MolToSmiles(final) == expected
           
    shutil.rmtree(est_file)


def test_ran_pol(rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder13')
    os.mkdir(est_file)
    os.chdir(est_file)
    
    mol_2=Chem.MolFromSmiles('BrCOCBr')
    mol_4=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
    mol_5=Chem.MolFromSmiles('BrCCBr')

    dia = Rnp(mol_2, mol_4,"Br").random_ab_copolymer(10, 0.4,force_field="MMFF94",relax_iterations=10)
    Fmt(dia).xyz_print("dia.xyz")
  
    tri = Rnp(mol_2, mol_4,"Br").random_abc_copolymer(mol_2, 15, 0.4, 0.6, force_field="MMFF94",relax_iterations=10)
    Fmt(tri).xyz_print("tri.xyz")

    b=Fld().seek_files("xyz")

    assert len(b) == 2
    
    shutil.rmtree(est_file)
