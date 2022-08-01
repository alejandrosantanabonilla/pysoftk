from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *

mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

a=Sm(mol_1,mol_2,"Br").monomer()

Fmt(a).xyz_print("first_monomer.xyz")

