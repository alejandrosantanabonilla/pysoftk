from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer.calculators import *

# Molecule 1 and 2 example
mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

# Creating a linear polymer
a=Sm(mol_1,mol_2,"Br")
k=a.mon_to_poly()
new=Lp(k,"Br",5,shift=1.0).linear_polymer(150)
Fmt(new).xyz_print("test_1.xyz")

# xtb optimisation
Opt("test_1.xyz").xtb_ff("xtb")
