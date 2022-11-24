from rdkit import Chem
from rdkit.Chem import AllChem

mol_1=Chem.MolFromSmiles('BrCOCBr')
mol_2=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')

from pysoftk.topologies.diblock import *

di=Db(mol_1,mol_2,"Br").diblock_copolymer(5,7,"MMFF",10)

from pysoftk.format_printers.format_mol import *

Fmt(di).xyz_print("diblock.xyz")
