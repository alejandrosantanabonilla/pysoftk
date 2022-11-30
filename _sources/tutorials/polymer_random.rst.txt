.. include:: /include/links.rst

.. _polymer_random:

=================================================================
Creating random Polymers
=================================================================

PySoftK_ supports the construction of random linear copolymers which consist of two or three different monomers. Our approach uses the command **Rnp(mol_1,mol_2,atom).random_ab_copolymer(len,pA,iter_ff,FF)** to build a random copolymer of length len that contains two monomers, **mol_1** and **mol_2**. The number of **mol_1** and **mol_2** monomers in the resulting polymer are determined from **pA len** and **(1-pA)len**, respectively.

.. literalinclude:: scripts/ran.py
   :lines: 1-5

Then, we can create a polymer with a diblock random topology, by using the following commands:

.. literalinclude:: scripts/ran.py
   :lines: 7-9

The new molecular unit (stored in the variable **ran**) can be used to replicate and form a diblock random polymer with a given desired number of units (in total **5**). A probability of **0.4** to be merged is defined by the user, while **10** steps for geometry optimization has been requested. The structure can be printed in XYZ format by adding the following lines of code as can be seen in this snipet:

.. literalinclude:: scripts/ran.py
   :lines: 11-13

By using a common visualization program (such as VMD_), the built structure **random2.xyz** can be displayed and the result as presented above

.. figure:: images/random2.png
   :align: center
   :figclass: align-center

   **Figure** Random diblock polymer with 5 units. 

Similarly PySoftK_ can be used to generate random linear copolymers consisting of three different monomers via the command **Rnp(mol_1,mol_2,atom).random_abc_copolymer(mol_3, len, pA, pB, iter_ff, FF)**. In this case, the resulting polymer will contain len monomers, such that the number of **mol_1**, **mol_2** and **mol_3** monomer is **pA len**, **pB len** and **(1-pA-pB) len**, respectively

.. literalinclude:: scripts/ran.py
   :lines: 17-21

The result of this function is displayed in the figure below.

.. figure:: images/random3.png
   :align: center
   :figclass: align-center

   **Figure** Random triblock 5-unit repetition polymer
