from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *
from pysoftk.htp_tools.calculator_htp import *

import os
from pathlib import Path
import pytest

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_htp_gfn_creation(rootdir):
  # Creating a tmp folder
  base_dir = os.path.join(rootdir, 'tmp_folder4')
  os.mkdir(base_dir)
  os.chdir(base_dir)

  # Test HTP for linear polymer chains created on the fly
  mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  a=Sm(mol_1,mol_2,"Br")

  molecules=[]
  for i in range(1,10):
    k=a.mon_to_poly()
    molecules.append(Lp(k,"Br",i,shift=1.0).linear_polymer("MMFF94",150*i))

  for idx, values in enumerate(molecules):
    Fmt(values).xyz_print("mono_"+str(idx)+".xyz")

  Fld().file_to_dir("xyz")
  xtb_identifier = "xtb"
  
  try:
    # Pass the command name or path to the constructor
    htp_calculator = Htp(xtb_command=xtb_identifier)

    # --- Run calculations ---
    print("\n--- Running GFN-XTB ---")
    htp_calculator.htp_xtb_gfn(
            directory=base_dir,
            max_work=1,
            num_cores=1,
            threshold="normal")


  except FileNotFoundError as e:
     print(f"\nInitialization Error: {e}")
     print("Please ensure 'xtb' is installed and in your PATH, or provide the full path.")
  except Exception as e:
     print(f"\nAn unexpected error occurred: {e}")
    
  working_dir=Path(os.getcwd())
  relaxed_str=[path for path in working_dir.glob("**/xtbopt.xyz")]

  p=Path(os.getcwd())
  a=[f for f in p.iterdir() if f.is_dir()]

  assert len(relaxed_str)==9
  assert len(a)==9

  shutil.rmtree(base_dir)

def test_htp_gff_creation(rootdir):
  # Creating a tmp folder
  base_dir = os.path.join(rootdir, 'tmp_folder70')
  os.mkdir(base_dir)
  os.chdir(base_dir)

  # Test HTP for linear polymer chains created on the fly
  mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  a=Sm(mol_1,mol_2,"Br")

  molecules=[]
  for i in range(1,10):
    k=a.mon_to_poly()
    molecules.append(Lp(k,"Br",i,shift=1.0).linear_polymer("MMFF94",150*i))

  for idx, values in enumerate(molecules):
    Fmt(values).xyz_print("mono_"+str(idx)+".xyz")

  Fld().file_to_dir("xyz")
  xtb_identifier = "xtb"
  
  try:
    # Pass the command name or path to the constructor
    htp_calculator = Htp(xtb_command=xtb_identifier)

    # --- Run calculations ---
    print("\n--- Running GFN-XTB ---")
    htp_calculator.htp_xtb_ff(
            directory=base_dir,
            max_work=1,
            num_cores=1,
            threshold="normal")


  except FileNotFoundError as e:
     print(f"\nInitialization Error: {e}")
     print("Please ensure 'xtb' is installed and in your PATH, or provide the full path.")
  except Exception as e:
     print(f"\nAn unexpected error occurred: {e}")
    
  working_dir=Path(os.getcwd())
  relaxed_str=[path for path in working_dir.glob("**/xtbopt.xyz")]

  p=Path(os.getcwd())
  a=[f for f in p.iterdir() if f.is_dir()]

  assert len(relaxed_str)==9
  assert len(a)==9

  shutil.rmtree(base_dir)


