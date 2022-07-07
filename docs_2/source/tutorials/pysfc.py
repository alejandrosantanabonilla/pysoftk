from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.linear_polymer.calculators import *
from pysoftk.format_printers.format_mol import *


# Molecule 1 and 2 example
mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

# Creating linear polymer
a=Sm(mol_1,mol_2,"Br")
k=a.mon_to_poly()
new=Lp(k,"Br",4,shift=1.25).linear_polymer(250)
Fmt(new).xyz_print("test_2.xyz")

# pyscf semiempirical geometry optimization
# "test_1.xyz" is provided as first argument in Opt
# the number of SCF cycles is set in pyscf_semi
Opt("test_2.xyz").pyscf_semi(1000)
