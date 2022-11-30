.. include:: /include/links.rst

.. _polymer_cyclic:

=================================================================
Creating cyclic Polymers
=================================================================

Constructing a polymer with cyclic topology can be carried out in PySoftK_ by employing the :mod:`pysoftk.topologies.ring` function. First, we need to import the initial molecule and read it using RDKit_, as shown in the following snippet:

.. literalinclude:: scripts/cyclic.py
   :lines: 1-5

Then, we can create a cyclic topology, by using the following commands:

.. literalinclude:: scripts/cyclic.py
   :lines: 7-9

The new molecular unit (stored in the variable **cyc**) can be used to replicate and form a polymer with a given desired number of units (in this case 8).

It is important to notice that we have used a **different** Force Field (Universal Force Field) and **increased** the number of iterations in the geometry optiomization process (iters=1000). This ensures a correct initial configuration highligthing the importance of always testing the correct parameters to generate correct structures.

To print the structure in XYZ format, one needs to add the following lines of conde as can be seen in this snipet:

.. literalinclude:: scripts/cyclic.py
   :lines: 11-13

By using a common visualization program (such as VMD_), the built structure **ring.xyz** can be displayed and the result as presented above

.. figure:: images/ring_2.png
   :width: 250px
   :align: center
   :height: 250px
   :figclass: align-center

   **Figure** Cyclic polymer with 8 repetition units. 
