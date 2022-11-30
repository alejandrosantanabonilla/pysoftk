.. include:: /include/links.rst

.. _polymer_diblock:

=================================================================
Creating diblock Polymers
=================================================================

For building linear diblock copolymers, PySoftK_ has a specific module (:mod:`pysoftk.topologies.diblock`)
which can be used. The command **Db(ma,mb,atom).diblock_copolymer(len_block_A,len_block_B,FF,iter_ff)**
is used to generate a diblock copolymer that consists of a block containing **len_block_A** monomers of **ma** and a block containing **len_block_B** monomers of **mb**. 


.. literalinclude:: scripts/diblock.py
   :lines: 1-5

Then, we can create a polymer with a diblock topology, by using the following commands:

.. literalinclude:: scripts/diblock.py
   :lines: 7-9

The new molecular unit (stored in the variable **di**) can be used to replicate and form a polymer with a given desired number of units (in total 12). The structure can be printed in XYZ format by adding the following lines of code as can be seen in this snipet:

.. literalinclude:: scripts/diblock.py
   :lines: 11-13

By using a common visualization program (such as VMD_), the built structure **diblock.xyz** can be displayed and the result as presented above

.. figure:: images/diblock.png
   :width: 350px
   :align: center
   :height: 350px
   :figclass: align-center

   **Figure** Diblock polymer with 12 units. 
