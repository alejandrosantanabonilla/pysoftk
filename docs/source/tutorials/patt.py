import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

mols=[Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'),
      Chem.MolFromSmiles('c1cc(oc1Br)Br'),
      Chem.MolFromSmiles('c1cc(sc1Br)Br')]

from pysoftk.topologies.diblock import *

string="ABC"
patt=Pt(string, mols, "Br").pattern_block_poly()

from pysoftk.format_printers.format_mol import *

Fmt(patt).xyz_print("patterned.xyz")

from itertools import permutations

alph1=permutations(["A", "B", "C"], 3)
string=[''.join(values) for idx, values in enumerate(alph1)]

for idx, values in enumerate(string):
    patt=Pt(values, mols, "Br").pattern_block_poly()
    Fmt(patt).xyz_print(values+".xyz")

from pysoftk.folder_manager.folder_creator import *
