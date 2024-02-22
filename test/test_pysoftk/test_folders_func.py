from pysoftk.folder_manager.folder_creator import *

import os
from pathlib import Path
import pathlib
import shutil
import pytest

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_seek_folders(rootdir):
    # Test 1:
    # Seek files with a given extension like: .smi, .mol, .pbd, etc
    # and print a list of the files.
    
    est_file = os.path.join(rootdir, 'tmp_folder1')
    os.mkdir(est_file)
    os.chdir(est_file)

    with open("mol_1.smi", "w") as f:
        f.write('c1cc(sc1Br)Br')

    with open("mol_2.smi", "w") as g:
        g.write('c1(ccc(cc1)Br)Br')

    a=Fld().seek_files("smi")
    assert len(a) == 2

    shutil.rmtree(est_file)

def test_organise_files(rootdir):
    # Test 2: Seek for .smi files and create the needed directories
    # using 2 cores.

    est_file = os.path.join(rootdir, 'tmp_folder2')
    os.mkdir(est_file)
    os.chdir(est_file)

    with open("mol_1.smi", "w") as f:
        f.write('c1cc(sc1Br)Br')

    with open("mol_2.smi", "w") as g:
        g.write('c1(ccc(cc1)Br)Br')
    
    Fld().file_to_dir("smi", 2)

    p=Path(os.getcwd())
    a=[f for f in p.iterdir() if f.is_dir()]
   
    assert len(a) == 2 

    shutil.rmtree(est_file)

def test_create_empty_dir(rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder3')
    os.mkdir(est_file)
    os.chdir(est_file)

    # Test 3: Create 2 folders with random names using 1 core.
    Fld().create(5)

    p=Path(os.getcwd())

    a=[f for f in p.iterdir() if f.is_dir()]
   
    assert len(a) == 5 

    shutil.rmtree(est_file)
