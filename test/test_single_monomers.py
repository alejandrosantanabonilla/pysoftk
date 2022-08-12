from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

from pysoftk.linear_polymer.super_monomer import *
from pysoftk.linear_polymer.linear_polymer import *

from pysoftk.format_printers.format_mol import *
from pysoftk.folder_manager.folder_creator import *

# SMILES all asymmetrical polymers

a=Chem.MolFromSmiles('c1cc(sc1CBr)Br')
#a=Chem.MolFromSmiles('BrOCCBr')
#a=Chem.MolFromSmiles('[C@@H](CBr)(Br)C(=O)N')
#a=Chem.MolFromSmiles('C([C@@H](CN)Br)Br')
#a=Chem.MolFromSmiles('C(/C=C(\C)/CBr)Br')
#a=Chem.MolFromSmiles('[C@@H](C(=O)OBr)(C)Br')
#a=Chem.MolFromSmiles('[C@@H](CBr)(Br)C(=O)OC')
#a=Chem.MolFromSmiles('[C@](CBr)(Br)(C(=O)OC)C')
#a=Chem.MolFromSmiles('c1c([C@@H](CBr)Br)cccc1')

#Very Problematic
#a=Chem.MolFromSmiles('[C@H](CBr)(OBr)C')
#b=Chem.MolToSmiles(a, isomericSmiles=False)
#c=Chem.MolFromSmiles(b)

# Original Embedding
AllChem.EmbedMolecule(a)

# New Embedding
#ps = rdDistGeom.ETKDGv3()
#AllChem.EmbedMolecule(c,ps)

new=Lp(a,"Br",7,shift=1.25).linear_polymer("MMFF",350)
Fmt(new).xyz_print("polymer_1_2.xyz")
