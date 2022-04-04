from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem import TorsionFingerprints


from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.mol_conformer import *
from pysoftk.format_printers.format_mol import *

# First set of Molecules for test 1
mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
mol_2=Chem.MolFromSmiles('c1cc(sc1Br)Br')

# Second set of molecules for test 2 
#mol_1=Chem.MolFromSmiles('c1cc(sc1Br)Br')
#mol_2=Chem.MolFromSmiles('c1(ccc(cc1)Br)Br')


# Original molecule
a=Sm(mol_1,mol_2,"Br")
k=a.monomer()
Fmt(k).xyz_print("old.xyz")

# Computing the conformer with highest energy
l=Mcon(k,1000,e_max=True).conformer()
Fmt(l).xyz_print("new.xyz")
