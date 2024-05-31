import pytest
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

from pysoftk.topologies.diblock import *
from pysoftk.topologies.ring import *
from pysoftk.topologies.branched import *
from pysoftk.topologies.ranpol import *

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *

from itertools import permutations
from pysoftk.tools.utils_rdkit import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

import os

# PYTESTS
testdata1=[('ABC', 'c1(ccccc1)c1ccc(o1)c1cccs1'),
           ('ACB', 'c1(ccccc1)c1ccc(s1)c1ccco1'),
           ('BAC', 'c1cc(oc1)c1ccc(cc1)c1cccs1'),
           ('BCA', 'c1cc(oc1)c1ccc(s1)c1ccccc1'),
           ('CAB', 'c1cc(sc1)c1ccc(cc1)c1ccco1'),
           ('CBA', 'c1cc(sc1)c1ccc(o1)c1ccccc1')]


testdata2=[('ABCD', 'c1cc(oc1)c1ccc(s1)c1ccc(cc1)CC'),
            ('ABDC', 'c1cc(oc1)c1ccc(s1)CCc1ccccc1'),
            ('ACBD', 'c1cc(oc1)c1ccc(cc1)c1ccc(s1)CC'),
            ('ACDB', 'c1cc(oc1)c1ccc(cc1)CCc1cccs1'),
            ('ADBC', 'c1cc(oc1)CCc1ccc(s1)c1ccccc1'),
            ('ADCB', 'c1cc(oc1)CCc1ccc(cc1)c1cccs1'),
            ('BACD', 'c1cc(sc1)c1ccc(o1)c1ccc(cc1)CC'),
            ('BADC', 'c1cc(sc1)c1ccc(o1)CCc1ccccc1'),
            ('BCAD', 'c1cc(sc1)c1ccc(cc1)c1ccc(o1)CC'),
            ('BCDA', 'c1cc(sc1)c1ccc(cc1)CCc1ccco1'),
            ('BDAC', 'c1cc(sc1)CCc1ccc(o1)c1ccccc1'),
            ('BDCA', 'c1cc(sc1)CCc1ccc(cc1)c1ccco1'),
            ('CABD', 'c1(ccccc1)c1ccc(o1)c1ccc(s1)CC'),
            ('CADB', 'c1(ccccc1)c1ccc(o1)CCc1cccs1'),
            ('CBAD', 'c1(ccccc1)c1ccc(s1)c1ccc(o1)CC'),
            ('CBDA', 'c1(ccccc1)c1ccc(s1)CCc1ccco1'),
            ('CDAB', 'c1(ccccc1)CCc1ccc(o1)c1cccs1'),
            ('CDBA', 'c1(ccccc1)CCc1ccc(s1)c1ccco1'),
            ('DABC', 'CCc1ccc(o1)c1ccc(s1)c1ccccc1'),
            ('DACB', 'CCc1ccc(o1)c1ccc(cc1)c1cccs1'),
            ('DBAC', 'CCc1ccc(s1)c1ccc(o1)c1ccccc1'),
            ('DBCA', 'CCc1ccc(s1)c1ccc(cc1)c1ccco1'),
            ('DCAB', 'CCc1ccc(cc1)c1ccc(o1)c1cccs1'),
            ('DCBA', 'CCc1ccc(cc1)c1ccc(s1)c1ccco1')]

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

@pytest.mark.parametrize("comb,expected", testdata1)
def test_threemols(comb, expected):
   """
   First test, 3 patterns
   """

   mols=[Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'),
         Chem.MolFromSmiles('c1cc(oc1Br)Br'),
         Chem.MolFromSmiles('c1cc(sc1Br)Br')]
   
   b=Pt(comb, mols, "Br").pattern_block_poly(swap_H=True)
   res=Chem.MolFromSmiles(b.write("smi"))

   assert Chem.MolToSmiles(res) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))

@pytest.mark.parametrize("comb,expected", testdata2)
def test_fourmols(comb, expected):
   """
   Second test, 4 patterns
   """
   mols=[Chem.MolFromSmiles('c1cc(oc1Br)Br'),
         Chem.MolFromSmiles('c1cc(sc1Br)Br'),
         Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'),
         Chem.MolFromSmiles('BrCCBr')]

   a=Pt(comb, mols, "Br").pattern_block_poly(swap_H=True)
   res=Chem.MolFromSmiles(a.write("smi"))
   assert Chem.MolToSmiles(res) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))

def test_bigamino(rootdir):
   """
   Third test, long polypeptide
   """
   #Loading precursors
   mol1tmp=os.path.join(rootdir, 'data/peg.mol')
   mol1=Chem.MolFromMolFile(mol1tmp, removeHs=False)

   mol2tmp=os.path.join(rootdir, 'data/pgla2.mol')
   mol2=Chem.MolFromMolFile(mol2tmp, removeHs=False)

   mols=[mol1, mol2]

   l1 = ["".join(["A" for i in range(45)])]
   l2 = ["".join(["B" for i in range(15)])]

   # Combine the lists and flatten the resulting list of lists
   pattern = "".join([item for sublist in [l1, l2] for item in sublist])

   model=Pt(str(pattern), mols, "Br").pattern_block_poly(relax_iterations=1500, force_field="MMFF", swap_H=True, rot_steps=1)
   res=Chem.MolFromSmiles(model.write("smi"))
   
   #Loading correct config
   #test_file = os.path.join(rootdir, 'data/test_lon_pp.mol')
   #final=Chem.MolFromMolFile(test_file)
   
   total_num=res.GetNumAtoms()
   assert int(total_num) == 285 

   

