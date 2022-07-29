from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *

import os
from pathlib import Path

def test_format_1():
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).xyz_print("test_2.xyz")
   b=Fld().seek_files("xyz")

   assert len(b) == 1
   os.remove("test_2.xyz")
   
def test_format_2():
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).pdb_print("test_2.mol")

   b=Fld().seek_files("pdb")
   assert len(b) == 1
   
   os.remove("test_2.pdb")
    
def test_format_3():
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).mol_print("test_2.mol")

   b=Fld().seek_files("mol")

   assert len(b) == 1
   
   os.remove("test_2.mol")
    
def test_sm_2():
   
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   
   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).pdb_print("test_2.pdb")

   b=Fld().seek_files("pdb")

   assert len(b) == 1
   
   os.remove("test_2.pdb")
   
def test_sm_3():

   mol_1=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   
   Fmt(a).mol_print("test_2.mol")
   b=Fld().seek_files("mol")

   assert len(b) == 1
   
   os.remove("test_2.mol")

def test_sm_4():
   mol_1=Chem.MolFromSmiles('c1c2c(ccc1Br)c1c(c3c2nc(cn3)Br)cccc1')
   mol_2=Chem.MolFromSmiles('n1nc(sc1Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()

   Fmt(a).mol_print("test_3.mol")
   b=Fld().seek_files("mol")

   assert len(b) == 1
   
   os.remove("test_3.mol")
