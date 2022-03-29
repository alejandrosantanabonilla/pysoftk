from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.format_printers.format_mol import *

# Molecules 1, hydrogen and embeding the molecule with rdkit
mol_1=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
m2=Chem.AddHs(mol_1)
AllChem.EmbedMolecule(m2)


# Print a smiles molecule to xyz format
Fmt(m2).xyz_print("m2.xyz")

# Print a smiles molecule to pdb format
#Fmt(m2).pdb_print("m2.pdb")

# Print a smiles molecule in mol format
#Fmt(m2).mol_print("m2.mol")



