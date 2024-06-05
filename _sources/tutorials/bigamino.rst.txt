.. include:: /include/links.rst

.. _bigamino:

=========================================================
Modelling soft-matter long chain polymers
=========================================================

The PySoftK_ function :mod:`pysoftk.topologies.diblock` has been developed to address the challenge of creating patterned polymers. In this new version of PySoftK_, we have enhanced the capabilities of this function for modelling long polymers. This tutorial explains how to use this function 


Modelling long chains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PySoftK_ enables the use of xtb_ as a **calculator** employing a customized Force Field implementation (called **xtb_ff**).
To introduce this functionality, one needs to import the corresponding RDKit_ auxiliary functions (line 1-2) and the needed modules from PySoftK_ as shown in the following snipet:

.. literalinclude:: scripts/large_amino.py
   :lines: 1-6

To demonstrate the functionality of this module, two aminoacids (PEG and PGLA) are used to create a *patterned-polymer* (and stored in the variable **mols**). The pattern can be created in python-lists containing user-defined sequences of letters. In this case we have used 45 *A's* and 15 *B's*. This is displayed in the following snipet: 

.. literalinclude:: scripts/large_amino.py
   :lines: 7-12

The script combines these lists and flattens them to create a single string (pattern) representing the entire block copolymer sequence as shown below:

.. literalinclude:: scripts/large_amino.py
   :lines: 13-15

The **Pt** function from PySoftK takes the pattern string (pattern), loaded molecules (mols), a placeholder atom ("Br"). The key parameters of this function includes
the relax_iterations (Number of relaxation steps for energy minimization), force_field (Force field used for energy calculations), swap_H
(allow hydrogen swapping during block placement) and rot_steps (number of conformers during placement). This is observed in the following snippet:
	   
.. literalinclude:: scripts/pysfc.py
   :lines: 15-17
	   
As shown in the previous example, the result of this operation is then printed in a **PDB** file called **test.pdb**, using the function **FMT(object)pdb_print(name_to_be_printed)**.    

.. literalinclude:: scripts/pysfc.py
   :lines: 17-19

The results of the geometrical optimization are printed to the file **test.pdb** and can be visualized using VMD_.

.. figure:: images/bigamino.png
   :width: 388px
   :align: center
   :height: 208px
   :figclass: align-center
   
