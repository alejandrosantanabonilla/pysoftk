from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.folder_manager.folder_creator import *
from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer.calculators import *

# Molecule 1 and 2 example
mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

# Creating linear polymer
a=Sm(mol_1,mol_2,"Br")
k=a.mon_to_poly()
new=Lp(k,"Br",4,shift=None).linear_polymer()
Fmt(new).xyz_print("mol_2.xyz")

# xtb optimisation
Opt("mol_2.xyz").xtb_ff("xtb")
