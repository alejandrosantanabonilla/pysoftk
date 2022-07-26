from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.mol_conformer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *
from pysoftk.torsional.torsional import *
from pysoftk.htp_tools.calculator_htp import *

import os
from pathlib import Path


def test_htp_gfn_creation():
=
def test_htp_ff_creation():
>>>>>>> 402967b925fa5e49c3ce90593c4d7f6a51e5a347

  # Test HTP for linear polymer chains created on the fly

  mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  a=Sm(mol_1,mol_2,"Br")

  molecules=[]
  for i in range(1,10):
    k=a.mon_to_poly()
<<<<<<< HEAD
    molecules.append(Lp(k,"Br",i,shift=1.0).linear_polymer("MMFF",150*i))
=======
    molecules.append(Lp(k,"Br",i,shift=1.0).linear_polymer(150*i))
>>>>>>> 402967b925fa5e49c3ce90593c4d7f6a51e5a347

  for idx, values in enumerate(molecules):
    Fmt(values).xyz_print("mono_"+str(idx)+".xyz")

  Fld().file_to_dir("xyz")

  # High-throughput calculations at the gfn-ff level of theory
<<<<<<< HEAD
  Htp("xyz").htp_xtb_gfn("xtb",4,1)
=======
  Htp("xyz").htp_xtb_ff("xtb",4,1)
>>>>>>> 402967b925fa5e49c3ce90593c4d7f6a51e5a347

  working_dir=Path(os.getcwd())
  relaxed_str=[path for path in working_dir.glob("**/xtbopt.xyz")]

  assert len(relaxed_str)==9

  p=Path(os.getcwd())
  a=[f for f in p.iterdir() if f.is_dir()]

  for i in a:
    shutil.rmtree(i)

<<<<<<< HEAD
    
=======
=======
def test_htp_gfn_creation():
>>>>>>> 402967b925fa5e49c3ce90593c4d7f6a51e5a347:test/test_htp_gfn.py


<<<<<<< HEAD:test/test_htp_ff.py
    
=======
  Fld().file_to_dir("xyz")

  # High-throughput calculations at the gfn-ff level of theory
  Htp("xyz").htp_xtb_gfn("xtb",4,1)

  working_dir=Path(os.getcwd())
  relaxed_str=[path for path in working_dir.glob("**/xtbopt.xyz")]

  assert len(relaxed_str)==9

  p=Path(os.getcwd())
  a=[f for f in p.iterdir() if f.is_dir()]

  for i in a:
    shutil.rmtree(i)
>>>>>>> 402967b925fa5e49c3ce90593c4d7f6a51e5a347:test/test_htp_gfn.py
>>>>>>> 402967b925fa5e49c3ce90593c4d7f6a51e5a347
