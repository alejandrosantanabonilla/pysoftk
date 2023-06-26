from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer.calculators import *

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *


def test_pyscf():

    # Molecule 1 and 2 example
    mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
    mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

    # Creating linear polymer
    a=Sm(mol_1,mol_2,"Br")
    k=a.mon_to_poly()
    new=Lp(k,"Br",2,shift=1.25).linear_polymer("MMFF94",250)
    Fmt(new).xyz_print("test_1.xyz")

    # pyscf semiempirical
    Opt("test_1.xyz").pyscf_semi(1000)

    a=Fld().seek_files("xyz")
  
    assert len(a) == 2

    os.remove("test_1.xyz")
    os.remove("pyscf_final.xyz")
