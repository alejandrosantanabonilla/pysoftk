from rdkit import Chem
from rdkit.Chem import AllChem

mol=Chem.MolFromSmiles("c1cc(ccc1Br)Br")
AllChem.EmbedMolecule(mol)

from pysoftk.topologies.ring import *

cyc=Rn(mol,'Br').pol_ring(8,FF="UFF",iters=1000)

from pysoftk.format_printers.format_mol import *

Fmt(cyc).xyz_print("ring.xyz")
