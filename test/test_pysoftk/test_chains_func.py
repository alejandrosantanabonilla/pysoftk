from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer.linear_polymer import Lpr

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *
from pysoftk.torsional.torsional import *
from pysoftk.torsional.mol_conformer import *

from pysoftk.folder_manager.folder_creator import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

import os
import pytest

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_lp_1(rootdir):
  #Sample molecules 1
  mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
  mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

  #Creating polymer 1
  a=Sm(mol_1,mol_2,"Br")
  k=a.mon_to_poly()

  new=Lp(k,"Br",15,shift=1.0).linear_polymer("MMFF94", relax_iterations=75)
  res_new=Chem.MolFromSmiles(new.write("smi"))

  mol1tmp=os.path.join(rootdir, 'data/polymer_1_2.mol')
  comp=Chem.MolFromMolFile(mol1tmp)
  
  assert Chem.MolToSmiles(comp) == Chem.MolToSmiles(res_new)

def test_lp_2(rootdir):

   #Sample molecule 2
   mol_3=Chem.MolFromSmiles('c1cc(sc1Br)Br')
   mol_4=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

   #Creating polymer 1
   r=Sm(mol_3,mol_4,"Br")
   g=r.mon_to_poly()

   new_1=Lp(g,"Br",15,shift=1.0).linear_polymer("MMFF94",relax_iterations=150)
   res_new=Chem.MolFromSmiles(new_1.write("smi"))   

   mol2tmp=os.path.join(rootdir, 'data/polymer_3_4.mol')
   comp=Chem.MolFromMolFile(mol2tmp)
   
   assert Chem.MolToSmiles(comp) == Chem.MolToSmiles(res_new)

   
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


def test_mol_conformer_1():
  conformer_generator = ConformerGenerator(num_conformers=5, forcefield="uff", make3D_steps=100, convergence=10)

  # Example usage 3: Overriding default parameters in ga_generate_conformers
  smiles_override = "CCO"  # SMILES string for ethanol
  base_name_override = "ethanol_ga_override"
  output_directory_override = "ethanol_ga_override_conformers"
  propanol_conformers_override = conformer_generator.ga_generate_conformers(
      smiles=smiles_override,
      num_conformers=8,       # Override default num_conformers
      forcefield="mmff94",    # Override default forcefield
      make3D_steps=50,        # Override default make3D_steps
      convergence=5           # Override default convergence
  )
  if propanol_conformers_override:
      conformer_generator.save_conformers_to_separate_files(propanol_conformers_override, base_name_override, output_dir=output_directory_override, file_format="xyz")
      print(f"Generated and saved {propanol_conformers_override.OBMol.NumConformers()} conformers (GA-like, overridden) in '{output_directory_override}'.")
  else:
      print("Failed to generate conformers for ethanol using GA-like method (overridden defaults).")


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



