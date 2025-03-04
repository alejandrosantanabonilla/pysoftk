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
  est_file = os.path.join(rootdir, 'tmp_folder4')
  os.mkdir(est_file)
  os.chdir(est_file)

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

  # High-throughput calculations at the gfn-ff level of theory
  Htp("xyz").htp_xtb_gfn("xtb",4,1)

  working_dir=Path(os.getcwd())
  relaxed_str=[path for path in working_dir.glob("**/xtbopt.xyz")]

  p=Path(os.getcwd())
  a=[f for f in p.iterdir() if f.is_dir()]

  assert len(relaxed_str)==9
  assert len(a)==9

  shutil.rmtree(est_file)

