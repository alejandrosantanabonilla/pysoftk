from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.format_printers.format_mol import *

# Molecules 1 and 2
mol_1=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')
mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')
<<<<<<< HEAD

# Creating a monomer with pysiftk
a=Sm(mol_1,mol_2,"Br").monomer()

# Print a monomer to xyz format
=======

# Creating a monomer with pysoftk
a=Sm(mol_1,mol_2,"Br").monomer()

# Print a smiles molecule to xyz format
>>>>>>> 6a6dcd42bb13973dd0e41db3bd1664581fb5429e
Fmt(a).xyz_print("test_2.xyz")

# Print a smiles molecule to pdb format
Fmt(a).pdb_print("test_2.pdb")

# Print a smiles molecule in mol format
Fmt(a).mol_print("test_2.mol")



