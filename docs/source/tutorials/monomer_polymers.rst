.. include:: /include/links.rst

.. _monomer_polymers:

=========================================================
From single molecules to polymers using PySoftK
=========================================================

PySoftK_ enables to create Polymers from a single monomer. To do so, the user needs to highligthing the positions where the merging process will take 
place, using a place-holder atom (as presented in `preparing_monomers`_).   

An example is in this section presented using a Styrene molecule. By using a common visualization program (such as VMD_), the functionalized molecule is 
presented below:

.. figure:: styrene.png
   :align: center
   :figclass: align-center

**(a)** Represents a functionalized single Styrene molecule in which two Bromine atoms are used as place-holders for bonding formation. **(b)** A Styrene 
polymer formed by 3 initial Styrene units. 

The process to build these kind of polymers is presented in this snapshot:

.. code-block:: python

   from rdkit import Chem
   from rdkit.Chem import AllChem
   from rdkit.Chem import rdDistGeom

   from pysoftk.linear_polymer.linear_polymer import *
   from pysoftk.format_printers.format_mol import *

   a=Chem.MolFromSmiles('c1c([C@@H](CBr)Br)cccc1')
   
   #Embedding is needed for being parsed as a pysoftk.object
   AllChem.EmbedMolecule(a)
   
   new=Lp(a,"Br",2,shift=1.25).linear_polymer("MMFF",350)
   Fmt(new).xyz_print("styrene_pol.xyz")
   
   
The Styrene molecule (**a**) is initially declared using SMILES format. The molecule has been embedded using one the methods available in RDKit_ and then 
parsed to *linear_polymer* order to create initial structures of a desired polymer.
   
