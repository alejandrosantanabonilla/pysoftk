from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.mol_conformer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *
from pysoftk.torsional.torsional import *
from pysoftk.folder_manager.folder_creator import *


def test_lp_1():

  #Sample molecules 1
  mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

  #Creating polymer 1
  a=Sm(mol_1,mol_2,"Br")
  k=a.mon_to_poly()

  new=Lp(k,"Br",15,shift=1.0).linear_polymer("MMFF94")
  Fmt(new).xyz_print("polymer_1_2.xyz")

  a=Fld().seek_files("xyz")
  
  assert len(a) == 1

  os.remove("polymer_1_2.xyz")
  

def test_lp_2():

   #Sample molecule 2
   mol_3=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_4=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   #Creating polymer 1
   r=Sm(mol_3,mol_4,"Br")
   g=r.mon_to_poly()

   new_1=Lp(g,"Br",15,shift=1.0).linear_polymer("MMFF94",150)
   Fmt(new_1).xyz_print("polymer_3_4.xyz")

   a=Fld().seek_files("xyz")
  
   assert len(a) == 1

   os.remove("polymer_3_4.xyz")
   
def test_mol_conformer():

  # Test 3 creating a conformer with MMFF force field
  a=Chem.MolFromSmiles('c1cc(oc1c1ccccc1)c1sccc1')
  m=AllChem.AddHs(a)

  # Computing the conformer with highest energy
  Mcon(m,10).conformer("conformers")

  mols=Chem.SDMolSupplier('conformers.sdf')
  res=[m for m in mols]
  
  assert len(res) == 10

  os.remove("conformers.sdf")

def test_torsional_list():

  #Test counting the torsional angles in a planar polymer
  
  molecules=Chem.MolFromSmiles('s1c(ccc1)c1sc(cc1)c1cccs1')
  lst_angles=Torsional(molecules).seek_angles()
  
  assert len(lst_angles) == 8

def test_torsional_images():

  #Test counting the torsional angles in a planar polymer
  
  molecules=Chem.MolFromSmiles('s1c(ccc1)c1sc(cc1)c1cccs1')
  Torsional(molecules).plot_trs_ang("mol_1")
  Torsional(molecules).plot_trs_ang("mol_1")

  a=Fld().seek_files("png")
  assert len(a) == 8
  
  for i in a:
    os.remove(i)
