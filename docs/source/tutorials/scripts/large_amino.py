from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

from pysoftk.topologies.diblock import *
from pysoftk.format_printers.format_mol import *

mols=[Chem.MolFromMolFile('peg.mol',removeHs=False), Chem.MolFromMolFile("pgla2.mol",removeHs=False)]

l1 = ["".join(["A" for i in range(45)])]
l2 = ["".join(["B" for i in range(10)])]


pattern = "".join([item for sublist in [l1, l2] for item in sublist])

a=Pt(str(pattern), mols, "Br").pattern_block_poly(relax_iterations=1500, force_field="MMFF", swap_H=True, rot_steps=1)

Fmt(a).pdb_print("test.pdb")
