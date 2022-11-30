.. include:: /include/links.rst
	     
.. _polymer_pattern:

=================================================================
Creating patterned Polymers
=================================================================

A Polymer can be constructed as a given combination of monomer(-s). Sometimes, there is the need to create combination of these units, which can be achieved in PySoftK_ by employing the :mod:`pysoftk.topologies.diblock` function. First, we need to import the initial molecules (stored in a Python_ list) using RDKit_, as shown in the following snippet:

.. literalinclude:: scripts/patt.py
   :lines: 1-7

Then, we can create a patterned topology, by defining a string that contains
a **pattern** and then parsed to the corresponding function as shown below:

.. literalinclude:: scripts/patt.py
   :lines: 9-12

In this case, we have used 3 monomers and they have been described an alphabetic order using the letters **ABC**, representing the order in which the moieties
will be merged.

The new molecular unit (stored in the variable **patt**) can be used to form a
patternd polymer. The structure can be printed in XYZ format by adding the following lines of code:

.. literalinclude:: scripts/patt.py
   :lines: 14-16

By using a common visualization program (such as VMD_), the built structure **ring.xyz** can be displayed and the result as presented above

.. figure:: images/ABC.png
   :width: 388px
   :align: center
   :height: 208px
   :figclass: align-center

   **Figure** Patterned polymer with 3 units (**ABC**). 


Combinatorial creation of patterned Polymers
===========================================================

If the user is interested in creating libraries with different combinations for a given set of monomers, this can done by PySoftK_ with standard Python_ libraries such as *itertools*. As an example, we are going to create all possible **permutations** that can be generated a set of 3 different monomers. This displayed in the following snippet:

.. literalinclude:: scripts/patt.py
   :lines: 18-25

This will create 6 different permutations **ABC/ACB/BAC/BCA/CAB/CBA** and all initial geometries are written in XYZ format. In the following section, we will present a functionality to organize files into folders developed in  PySoftK_.

.. list-table:: 
   :class: borderless

   * - .. image:: images/ABC.png
     - .. image:: images/ACB.png 
   

.. list-table:: 
   :class: borderless
       
   * - .. image:: images/BAC.png
     - .. image:: images/BCA.png
   
   
.. list-table:: 
   :class: borderless

   * - .. image:: images/CAB.png
     - .. image:: images/CBA.png 

   
