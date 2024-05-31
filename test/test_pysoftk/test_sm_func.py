import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *
import os

testdata1=[('c1(ccc(cc1)Br)Br',
            '[H]c1sc(-c2c([H])c([H])c([H])c([H])c2[H])c([H])c1[H]')]

testdata2=[('c1cc(sc1Br)Br',
            '[H]c1sc(-c2sc([H])c([H])c2[H])c([H])c1[H]')]

testdata3=[('c1(ccc(cc1)Br)Br',
            'c1ccc(-c2ccccc2)cc1')]

testdata4=[('c1c2c(ccc1Br)c1c(c3c2nc(cn3)Br)cccc1',
            'c1cnc2c3ccccc3c3ccc(c4nncs4)cc3c2n1')]

testdata5=[('c1c(ccc2c1S(=O)(=O)c1c(C2=O)ccc(c1)Br)Br',
           'O=C1c2ccccc2S(=O)(=O)c2cc(-c3ccccc3)ccc21')]

@pytest.mark.parametrize("comb,expected", testdata1)
def test_sm_1(comb, expected):
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles(comb)

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)

   a=Chem.CanonSmiles(test_p)
   b=Chem.CanonSmiles(expected)

   assert a==b

@pytest.mark.parametrize("comb,expected", testdata2)
def test_sm_2(comb, expected):
   
   mol_1=Chem.MolFromSmiles(comb)
   mol_2=Chem.MolFromSmiles(comb)
   
   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)   
   
   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(expected)

   assert a==b   

@pytest.mark.parametrize("comb,expected", testdata3)
def test_sm_3(comb, expected):

   mol_1=Chem.MolFromSmiles(comb)
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)

   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(expected)

   assert a==b 

@pytest.mark.parametrize("comb,expected", testdata4)
def test_sm_4(comb, expected):
   
   mol_1=Chem.MolFromSmiles(comb)
   mol_2=Chem.MolFromSmiles('n1nc(sc1Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)
   
   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(expected)

   assert a==b

@pytest.mark.parametrize("comb,expected", testdata5)
def test_sm_5(comb, expected):
   mol_1=Chem.MolFromSmiles(comb)
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)

   #result=
   
   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(expected)

   assert a==b  



