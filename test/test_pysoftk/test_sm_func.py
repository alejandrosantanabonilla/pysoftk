from rdkit import Chem
from rdkit.Chem import AllChem
from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *

def test_sm_1():
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)
   result='[H]c1sc(-c2c([H])c([H])c([H])c([H])c2[H])c([H])c1[H]'

   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(result)

   assert a==b

def test_sm_2():
   
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   
   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)   
   result='[H]c1sc(-c2sc([H])c([H])c2[H])c([H])c1[H]'

   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(result)

   assert a==b   

def test_sm_3():

   mol_1=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)

   result = 'c1ccc(-c2ccccc2)cc1'

   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(result)

   assert a==b 

def test_sm_4():
   mol_1=Chem.MolFromSmiles('c1c2c(ccc1Br)c1c(c3c2nc(cn3)Br)cccc1')
   mol_2=Chem.MolFromSmiles('n1nc(sc1Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)

   result='c1cnc2c3ccccc3c3ccc(c4nncs4)cc3c2n1'
   
   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(result)

   assert a==b
   
def test_sm_5():
   mol_1=Chem.MolFromSmiles('c1c(ccc2c1S(=O)(=O)c1c(C2=O)ccc(c1)Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   temp=Sm(mol_1,mol_2,"Br").monomer()
   test_p=Chem.MolToSmiles(temp)

   result='O=C1c2ccccc2S(=O)(=O)c2cc(-c3ccccc3)ccc21'
   
   a = Chem.CanonSmiles(test_p)
   b = Chem.CanonSmiles(result)

   assert a==b  



