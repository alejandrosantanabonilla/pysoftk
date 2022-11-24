from rdkit import Chem
from rdkit.Chem import AllChem

core=Chem.MolFromSmiles("BrN(Br)CCN(Br)Br")
arm=Chem.MolFromSmiles("C(CC(=O)OCCOC(=O)CCBr)Br")

from pysoftk.topologies.branched import *

bran=Bd(core,arm,"Br").branched_polymer()

from pysoftk.format_printers.format_mol import *
Fmt(bran).xyz_print("branched.xyz")
