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

def test_file_to_dir(rootdir):
    """
    Tests the file_to_dir function.
    """
    test_dir = os.path.join(rootdir, 'tmp_folder200')
    os.mkdir(test_dir)
    os.chdir(test_dir)

    # Create test files
    with open("mol_1.smi", "w") as f:
        f.write('c1cc(sc1Br)Br')

    with open("mol_2.smi", "w") as g:
        g.write('c1(ccc(cc1)Br)Br')

    # Generate fixed names
    fixed_names = Fld().fxd_name("molecule", 2)  

    # Move files
    Fld().file_to_dir(format_extension="smi", num_cores=2, fixed_names=fixed_names)

    p=Path(os.getcwd())
    a=[f for f in p.iterdir() if f.is_dir()]

    assert len(a) == 2 

    # Cleanup
    shutil.rmtree(test_dir)


def test_create_empty_dir(rootdir):

    est_file = os.path.join(rootdir, 'tmp_folder3')
    os.mkdir(est_file)
    os.chdir(est_file)

    # Test 3: Create 5 folders with random names using 1 core.
    Fld().create(5)

    p=Path(os.getcwd())

    a=[f for f in p.iterdir() if f.is_dir()]
   
    assert len(a) == 5 

    shutil.rmtree(est_file)

def test_move_files_to_folder(rootdir):
    # New test using temp_directory fixture
    """Test the move_files_to_folder function."""

    temp_dir_path = os.path.join(rootdir, "my_temp_dir_1000")
    os.makedirs(temp_dir_path)
    os.chdir(temp_dir_path)
    
    # Create some test files
    with open("file1.txt", "w") as f:
        f.write("Test file 1")
    with open("file2.txt", "w") as f:
        f.write("Test file 2")
    with open("other_file.dat", "w") as f:  # A file with a different extension
        f.write("Not a text file")

    # Instantiate the Fld class
    fld = Fld()

    # Call the function to move the files
    fld.move_files_to_folder("txt", "text_files")

    # Assertions to check if the files were moved correctly
    text_files_dir = os.path.join(temp_dir_path, "text_files")

    # Use temp_directory
    assert os.path.exists(text_files_dir)
    assert os.path.exists(os.path.join(text_files_dir, "file1.txt"))
    assert os.path.exists(os.path.join(text_files_dir, "file2.txt"))
    assert os.path.exists(os.path.join(temp_dir_path, "other_file.dat"))

    shutil.rmtree(temp_dir_path)
