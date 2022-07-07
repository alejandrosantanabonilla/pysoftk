from pysoftk.torsional.torsional import *

from rdkit import Chem
from rdkit.Chem import AllChem

#### Torsional angle detection for a molecule in a *.smi file ###

molecules=Chem.MolFromSmiles('s1c(ccc1)c1sc(cc1)c1cccs1')

# This list records all tuples of atoms which are involved in
# torsional angles for linear polymers.

lst_angles=Torsional(molecules).seek_angles()


#lst_tor_angl=seek_pol_tor_angles(molecules)

# This for loop iterates over all the previous cases
# and creates pictures of the molecule and highlights
# the atoms of every detected torsional angle

Torsional(molecules).plot_trs_ang("mol_1")
