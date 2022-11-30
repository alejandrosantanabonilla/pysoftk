from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *

from pysoftk.htp_tools.calculator_htp import *

from rdkit import Chem
from rdkit.Chem import AllChem


# Test HTPc tools for linear polymer chains created on the fly

mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
a=Sm(mol_1,mol_2,"Br")

molecules=[]
for i in range(1,10):
   k=a.mon_to_poly()
   molecules.append(Lp(k,"Br",i,shift=1.0).linear_polymer(150*i))

for idx, values in enumerate(molecules):
   Fmt(values).xyz_print("mono_"+str(idx)+".xyz")

a=Fld().seek_files("xyz")
Fld().file_to_dir("xyz",2)

# High-throughput calculations at the gfn-ff level of theory
Htp("xyz").htp_xtb_ff("xtb",4,1)
