from pysoftk.folder_manager.folder_creator import *

# Create 2 folders with unique labels using 1 core.
Fld().create(2)

from rdkit import Chem
from rdkit.Chem import AllChem

print('c1cc(sc1Br)Br',file=open("mol_1.smi", 'w'))
print('c1(ccc(cc1)Br)Br', file=open("mol_2.smi",'w'))

# Seek files with a given extension like: .smi, .mol, .pbd, .xyz, etc
a=Fld().seek_files("smi")

from rdkit import Chem
from rdkit.Chem import AllChem

print('c1cc(sc1Br)Br',file=open("mol_1.smi", 'w'))
print('c1(ccc(cc1)Br)Br', file=open("mol_2.smi",'w'))

# Test 3: Seek for .smi files and create the needed directories using 2 cores.
Fld().file_to_dir("smi",2)


