from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.mol_conformer import *
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *

#Sample molecules 1
mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

#Creating polymer 1
a=Sm(mol_1,mol_2,"Br")
k=a.mon_to_poly()
new=Lp(k,"Br",3).linear_polymer(10)
Fmt(new).xyz_print("polymer_1_2.xyz")

#Sample molecule 2
#mol_3=Chem.MolFromSmiles('c1cc(sc1Br)Br')
#mol_4=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

#Creating polymer 1
#r=Sm(mol_3,mol_4,"Br")
#g=r.mon_to_poly()
#new_1=Lp(g,"Br",5).linear_polymer(150)
#Fmt(new_1).xyz_print("polymer_3_4.xyz")
