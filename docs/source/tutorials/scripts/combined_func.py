import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import permutations

from pysoftk.topologies.diblock import *
from pysoftk.linear_polymer.utils import *

# Creating permutations for 4 different units
alph1=permutations(["A", "B", "C", "D"], 4)
string=[''.join(values) for idx, values in enumerate(alph1)]

mols=[Chem.MolFromSmiles('c1cc(oc1Br)Br'),
      Chem.MolFromSmiles('c1cc(sc1Br)Br'),
      Chem.MolFromSmiles('c1(ccc(cc1)Br)Br'),
      Chem.MolFromSmiles('c1(ccc(nc1)Br)Br')]

patt=Pt(string[0], mols, "Br").pattern_block_poly(swap_H=False)

# Linear polymer module 
from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *

new=Lp(patt,"Br", 2, shift=1.0).linear_polymer("MMFF", 5)
Fmt(new).xyz_print("ABCD_pol.xyz")
