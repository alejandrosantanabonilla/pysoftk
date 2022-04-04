from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.format_printers.format_mol import *

# Molecules 1 and 2
mol_1=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')

# Print a smiles molecule to xyz format
Fmt(mol_1).xyz_print("test_2.xyz")

# Print a smiles molecule to pdb format
#Fmt(mol_1).pdb_print("test_2.pdb")

# Print a smiles molecule in mol format
#Fmt(mol_1).mol_print("test_2.mol")



