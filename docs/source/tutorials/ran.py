from rdkit import Chem
from rdkit.Chem import AllChem

mol_1=Chem.MolFromSmiles("BrCOCBr")
mol_2=Chem.MolFromSmiles("c1(ccc(cc1)Br)Br")

from pysoftk.topologies.ranpol import *

ran=Rnp(mol_1,mol_2,"Br").random_ab_copolymer(5,0.4,10)

from pysoftk.format_printers.format_mol import *

Fmt(ran).xyz_print("random2.xyz")

# triblock random polymer

from pysoftk.topologies.ranpol import *

mol_3=Chem.MolFromSmiles("c1cc(oc1Br)Br")
tri=Rnp(mol_1,mol_2,"Br").random_abc_copolymer(mol_3,5,0.4,0.2,10)
Fmt(tri).xyz_print("random3.xyz")


