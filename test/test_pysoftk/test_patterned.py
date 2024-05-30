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

# PYTESTS


testdata1=[('ABC', '[H]c1c(-c2sc(Br)c([H])c2[H])oc(-c2c([H])c([H])c(Br)c([H])c2[H])c1[H]'),
           ('ACB', '[H]c1c(Br)oc(-c2sc(-c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('BAC', '[H]c1c(Br)oc(-c2c([H])c([H])c(-c3sc(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('BCA', '[H]c1c(Br)oc(-c2sc(-c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('CAB', '[H]c1c(Br)oc(-c2c([H])c([H])c(-c3sc(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('CBA', '[H]c1c(-c2sc(Br)c([H])c2[H])oc(-c2c([H])c([H])c(Br)c([H])c2[H])c1[H]')]


testdata2=[('ABCD', '[H]c1c(Br)oc(-c2sc(-c3c([H])c([H])c(C([H])([H])C([H])([H])Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('ABDC', '[H]c1c(Br)oc(-c2sc(C([H])([H])C([H])([H])c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('ACBD', '[H]c1c(Br)oc(-c2c([H])c([H])c(-c3sc(C([H])([H])C([H])([H])Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('ACDB', '[H]c1c(Br)oc(-c2c([H])c([H])c(C([H])([H])C([H])([H])c3sc(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('ADBC', '[H]c1c(Br)oc(C([H])([H])C([H])([H])c2sc(-c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('ADCB', '[H]c1c(Br)oc(C([H])([H])C([H])([H])c2c([H])c([H])c(-c3sc(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('BACD', '[H]c1c(-c2sc(Br)c([H])c2[H])oc(-c2c([H])c([H])c(C([H])([H])C([H])([H])Br)c([H])c2[H])c1[H]'),
           ('BADC', '[H]c1c(-c2sc(Br)c([H])c2[H])oc(C([H])([H])C([H])([H])c2c([H])c([H])c(Br)c([H])c2[H])c1[H]'),
           ('BCAD', '[H]c1c(-c2c([H])c([H])c(-c3sc(Br)c([H])c3[H])c([H])c2[H])oc(C([H])([H])C([H])([H])Br)c1[H]'),
           ('BCDA', '[H]c1c(Br)oc(C([H])([H])C([H])([H])c2c([H])c([H])c(-c3sc(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('BDAC', '[H]c1c(-c2c([H])c([H])c(Br)c([H])c2[H])oc(C([H])([H])C([H])([H])c2sc(Br)c([H])c2[H])c1[H]'),
           ('BDCA', '[H]c1c(Br)oc(-c2c([H])c([H])c(C([H])([H])C([H])([H])c3sc(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('CABD', '[H]c1c(-c2sc(C([H])([H])C([H])([H])Br)c([H])c2[H])oc(-c2c([H])c([H])c(Br)c([H])c2[H])c1[H]'),
           ('CADB', '[H]c1c(-c2c([H])c([H])c(Br)c([H])c2[H])oc(C([H])([H])C([H])([H])c2sc(Br)c([H])c2[H])c1[H]'),
           ('CBAD', '[H]c1c(-c2sc(-c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])oc(C([H])([H])C([H])([H])Br)c1[H]'),
           ('CBDA', '[H]c1c(Br)oc(C([H])([H])C([H])([H])c2sc(-c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('CDAB', '[H]c1c(-c2sc(Br)c([H])c2[H])oc(C([H])([H])C([H])([H])c2c([H])c([H])c(Br)c([H])c2[H])c1[H]'),
           ('CDBA', '[H]c1c(Br)oc(-c2sc(C([H])([H])C([H])([H])c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('DABC', '[H]c1c(-c2sc(-c3c([H])c([H])c(Br)c([H])c3[H])c([H])c2[H])oc(C([H])([H])C([H])([H])Br)c1[H]'),
           ('DACB', '[H]c1c(-c2c([H])c([H])c(-c3sc(Br)c([H])c3[H])c([H])c2[H])oc(C([H])([H])C([H])([H])Br)c1[H]'),
           ('DBAC', '[H]c1c(-c2sc(C([H])([H])C([H])([H])Br)c([H])c2[H])oc(-c2c([H])c([H])c(Br)c([H])c2[H])c1[H]'),
           ('DBCA', '[H]c1c(Br)oc(-c2c([H])c([H])c(-c3sc(C([H])([H])C([H])([H])Br)c([H])c3[H])c([H])c2[H])c1[H]'),
           ('DCAB', '[H]c1c(-c2sc(Br)c([H])c2[H])oc(-c2c([H])c([H])c(C([H])([H])C([H])([H])Br)c([H])c2[H])c1[H]'),
           ('DCBA', '[H]c1c(Br)oc(-c2sc(-c3c([H])c([H])c(C([H])([H])C([H])([H])Br)c([H])c3[H])c([H])c2[H])c1[H]')]


@pytest.mark.parametrize("comb,expected", testdata1)
def test_threemols(comb, expected):
   """
   First test, 3 patterns
   """

   mols=[Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'),
         Chem.MolFromSmiles('c1cc(oc1Br)Br'),
         Chem.MolFromSmiles('c1cc(sc1Br)Br')]

   b=Pt(comb, mols, "Br").pattern_block_poly(swap_H=False)
   assert Chem.MolToSmiles(b) == expected


@pytest.mark.parametrize("comb,expected", testdata2)
def test_fourmols(comb, expected):
   """
   Second test, 4 patterns
   """
   mols=[Chem.MolFromSmiles('c1cc(oc1Br)Br'),
         Chem.MolFromSmiles('c1cc(sc1Br)Br'),
         Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'),
         Chem.MolFromSmiles('BrCCBr')]


   a=Pt(str(comb), mols, "Br").pattern_block_poly(swap_H=False)
   assert Chem.MolToSmiles(a) == expected





   

