from pysoftk.folder_manager.folder_creator import *
import os
from pathlib import Path

def test_seek_folders():

    # Test 1:
    # Seek files with a given extension like: .smi, .mol, .pbd, etc
    # and print a list of the files.
    
    print('c1cc(sc1Br)Br',file=open("mol_1.smi", 'w'))
    print('c1(ccc(cc1)Br)Br', file=open("mol_2.smi",'w'))
    
    a=Fld().seek_files("smi")

    assert len(a) == 2

    os.remove("mol_1.smi")
    os.remove("mol_2.smi")


def test_organise_files():

   print('c1cc(sc1Br)Br',file=open("mol_1.smi", 'w'))
   print('c1(ccc(cc1)Br)Br', file=open("mol_2.smi",'w'))
    
   # Test 2: Seek for .mol files and create the needed directories
   # using 2 cores.

   Fld().file_to_dir("smi",2)

   p=Path(os.getcwd())
   shutil.rmtree(os.path.join(p, "__pycache__"))

   a=[f for f in p.iterdir() if f.is_dir()]
   
   assert len(a) == 2 

   # Remove directories to perform next test
   for i in a:
       shutil.rmtree(i)


def test_crete_empty_dir():
    
  # Test 3: Create 2 folders with random names using 1 core.
  Fld().create(2,1)

  p=Path(os.getcwd())

  a=[f for f in p.iterdir() if f.is_dir()]
   
  assert len(a) == 2 

   # Remove directories to perform next test
  for i in a:
      shutil.rmtree(i)  
