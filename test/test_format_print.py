from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *
from pysoftk.linear_polymer.linear_polymer import *

import os
from pathlib import Path

import pytest

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_format_1(rootdir):
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()

   est_file = os.path.join(rootdir, 'tmp_folder5')
   os.mkdir(est_file)
   os.chdir(est_file)
   
   Fmt(a).xyz_print("test_2.xyz")
   b=Fld().seek_files("xyz")

   assert len(b) == 1
   shutil.rmtree(est_file)
   
def test_format_2(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder5')
   os.mkdir(est_file)
   os.chdir(est_file)
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).pdb_print("test_2.pdb")

   b=Fld().seek_files("pdb")
   assert len(b) == 1
   
   shutil.rmtree(est_file)
    
def test_format_3(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder6')
   os.mkdir(est_file)
   os.chdir(est_file)
    
   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).mol_print("test_2.mol")

   b=Fld().seek_files("mol")

   assert len(b) == 1
   
   shutil.rmtree(est_file)

def test_format_4(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder7')
   os.mkdir(est_file)
   os.chdir(est_file) 

   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

   AllChem.EmbedMolecule(mol_1)
   AllChem.EmbedMolecule(mol_2)

   #Creating polymer 1
   a=Sm(mol_1,mol_2,"Br")
   k=a.mon_to_poly()
   new=Lp(k,"Br", 3, shift=1.0).linear_polymer("MMFF94")

   #Different formats
   Fmt(new).format_print("molecule_3.mol2")
   Fmt(new).format_print("molecule_3.smi")
   Fmt(new).format_print("molecule_3.can")

   b=Fld().seek_files("mol2")
   assert len(b) == 1


   c=Fld().seek_files("smi")
   assert len(c) == 1
   

   d=Fld().seek_files("can")
   assert len(d) == 1

   shutil.rmtree(est_file)

def test_format_5(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder8')
   os.mkdir(est_file)
   os.chdir(est_file)

   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

   AllChem.EmbedMolecule(mol_1)
   AllChem.EmbedMolecule(mol_2)

   #Creating polymer 1
   a=Sm(mol_1,mol_2,"Br")
   k=a.mon_to_poly()
   new=Lp(k,"Br", 3, shift=1.0).linear_polymer("MMFF94")

   #Printing unit2
   Fmt(new).format_print("unit.mol")

   #Converting to mol2
   Cnv("unit.mol").file_in_out()

   #Checking the formats
   b=Fld().seek_files("mol2")
   assert len(b) == 1
   
   shutil.rmtree(est_file)
   
def test_sm_2(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder9')
   os.mkdir(est_file)
   os.chdir(est_file)

   mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   
   a=Sm(mol_1,mol_2,"Br").monomer()
   Fmt(a).pdb_print("test_2.pdb")

   b=Fld().seek_files("pdb")

   assert len(b) == 1
   
   shutil.rmtree(est_file)
   
def test_sm_3(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder10')
   os.mkdir(est_file)
   os.chdir(est_file)

   mol_1=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
   mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()
   
   Fmt(a).mol_print("test_2.mol")
   b=Fld().seek_files("mol")

   assert len(b) == 1
   
   shutil.rmtree(est_file)

def test_sm_4(rootdir):

   est_file = os.path.join(rootdir, 'tmp_folder11')
   os.mkdir(est_file)
   os.chdir(est_file)

   mol_1=Chem.MolFromSmiles('c1c2c(ccc1Br)c1c(c3c2nc(cn3)Br)cccc1')
   mol_2=Chem.MolFromSmiles('n1nc(sc1Br)Br')

   a=Sm(mol_1,mol_2,"Br").monomer()

   Fmt(a).mol_print("test_3.mol")
   b=Fld().seek_files("mol")

   assert len(b) == 1
   
   shutil.rmtree(est_file)

   
