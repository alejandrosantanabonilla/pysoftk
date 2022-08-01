from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *

from pysoftk.linear_polymer.mol_conformer import *

#Sample molecules 1
mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

a=Sm(mol_1,mol_2,"Br")
k=a.mon_to_poly()

new=Lp(k,"Br",3,shift=1.0).linear_polymer()
Fmt(new).xyz_print("polymer_1_2.xyz")
