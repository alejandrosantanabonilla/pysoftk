.. include:: /include/links.rst

.. _calculator:

=========================================================
Computing using **ab-initio** calculators
=========================================================

The PySoftK_ function :mod:`pysoftk.linear\_polymer.calculators` has been developed to perform geometry optimization for a single system or for multiple systems as a part of a high-throughput calculation workflow. 

The following projects have been added to PySoftK_ to perform geometry optimization calculations at different levels of theory.


Geometry optimization using GFN-FF 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PySoftK_ enables the use of xtb_ as a **calculator** employing a customized Force Field implementation (called **xtb_ff**).
To introduce this functionality, one needs to import the corresponding RDKit_ auxiliary functions (line 1-2) and the needed modules from PySoftK_ as shown in the following snipet:

.. literalinclude:: scripts/xtb_ff.py
   :lines: 1-7

To demonstrate the functionality of this module, two molecules (Benzene and Thiol) are used to create a *super-monomer* (and stored in the variable **a**). The object **a** is later used as an initial polymer unit (variable **k**) and replicated 5 times (and stored in variable **new**). The result of this operation is then printed in a **XYZ** file called **test_1.xyz**, using the function **FMT(object).xyz_print(name_to_be_printed)**.

.. literalinclude:: scripts/xtb_ff.py
   :lines: 9-17

The calculator module :mod:`pysoftk.linear\_polymer.calculators.Opt` can be used to perform a geometry optimization by providing a valid **XYZ** file 

.. warning::
   
   - The only xtb_ accepted geometry format is **XYZ**. 
   - To globally enable the xtb_ executable, one can define the **PATH** 
     and the **XTBHOME** variables using the commands:

     export PATH="$PATH:$HOME/PATH_TO_FOLDER/xtb-X.X.X/bin"
     export XTBHOME="PATH_TO_XTB_FOLDER/xtb-X.X.X"

   
The *Opt* module receives the name of the **XYZ** to be used by xtb_ whilst the function *xtb_ff* accepts the **absolute path/alias** where the xtb_ executable is stored (see note for creating an *alias* for xtb_ executable).

.. literalinclude:: scripts/xtb_ff.py
   :lines: 19-20

The results of the geometrical optimization are printed to the file **xtbopt.xyz** and can be visualized using VMD_.
	   
Geometry optimization using PySCF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PySoftK_ also permits the use of pyscf_ as a **calculator** at the semiempirical-mindo3 level of theory. Initially, one needs to import the corresponding PySoftK_ modules and RDKit_ auxiliary functions as shown below: 

.. literalinclude:: scripts/pysfc.py
   :lines: 1-7

As shown in the previous example, the result of this operation is then printed in a **XYZ** file called **test_2.xyz**, using the function **FMT(object).xyz_print(name_to_be_printed)**.    

.. literalinclude:: scripts/pysfc.py
   :lines: 10-18

.. note:: The only pyscf_ accepted format is **XYZ**.
	   
As displayed in the previous section, the *Opt* module accepts as an argument the name of the **xyz** file whilst the *pysc_semi*  receives the user-defined number of SCF cycles. 

.. literalinclude:: scripts/pysfc.py
   :lines: 20-23

The result of the geometrical optimization is printed to the file **pyscf_final.xyz** and can be visualized using VMD_.
	   
High-Throughput Calculations using the calculators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PySoftK_ has been designed with the aim of allowing high-throughput calculations (HTPc) at different levels of theory. This is achieved by using the PySoftK_ module :mod:`pysoftk.htp\_tools.calculator\_htp.Htp`, which can be combined with :mod:`pysoftk.folder\_manager.folder\_creator.Fld` for setting up and creating numerical workflows that include *ab-initio* calculations.

The needed PySoftK_ modules and RDKit_ functions are imported as displayed in the following lines of code:

.. literalinclude:: scripts/htp_test.py
   :lines: 1-10

For this example, a set of 10 linear Polymers (using Thiol as initial monomer) are assembled and merged into one **super-monomer**. 

.. literalinclude:: scripts/htp_test.py
   :lines: 19-22

Each of these new molecules is uniquely printed using the name **mono_IDX.xyz**
(employing the module :mod:`pysoftk.format\_printers package`) as indicated in the following lines of code:

.. literalinclude:: scripts/htp_test.py
   :lines: 24-25

The resulting files are organized using the module :mod:`pysoftk.folder\_manager.folder\_creator.Fld` where the function *seek_files* browses for files with the extension **.xyz**. To organize the files into folders, PySoftK_ utilizes the function *file_to_dir* indicating the type of files (**.xyz**) and the number of cores (**2**) which are requested to perform this operation:

.. literalinclude:: scripts/htp_test.py
   :lines: 27-28

The resulting directories can be browsed using the bash command **tree**

.. code-block:: bash

   $ tree
   .
   ├── 1b03a28c7fc94cb081071ad04947957a
   │   └── mono_4.xyz
   ├── 1ebe629956554d31a1ffcf141fdb6923
   │   └── mono_7.xyz
   ├── 233a3381070640088bc983ebd8c91d5c
   │   └── mono_0.xyz
   ├── 565413cfde9747548c07430c4d2ea2db
   │   └── mono_6.xyz
   ├── 62338a8a73c54d1daef85a19919cd4e2
   │   └── mono_3.xyz
   ├── 63f899b660504aa9b4ee4854712f876e
   │   └── mono_2.xyz
   ├── 67552ebcc87a496e85ce8f8198aa6af3
   │   └── mono_5.xyz
   ├── c6a52c74b4eb414f93837643d2d49736
   │   └── mono_1.xyz
   ├── e7518a8bf80540a6a3bc49567788db7f
   │   └── mono_8.xyz
   └── htp_test.py


The module :mod:`pysoftk.htp\_tools.calculator\_htp.Htp` enables the user to perform HTPc in a simple manner. Meanwhile, the command **HTP(format).htp_xtb_ff(executable, processes, cores)** can be used to perform geometry optimizations at different theory levels accessible in PySoftK_, as is indicated in the following lines of code:

.. literalinclude:: scripts/htp_test.py
   :lines: 30-31

PySoftK_ HTPc tools have been developed to perform parallel calculations based on the structures created in the folders. This is achieved by using asyncio_ parallel paradigm. The results of the workflow can be browsed using **tree** bash command as shown in the following snipet:

.. code-block:: bash

   $ tree
   .
   ├── 0d20dd8cab29474dbdcb7364c22cf3fb
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_4.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── 0e5a6600dbcd49e48a220ec104065bf0
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_7.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── 2da5bca994c84c318a1f46aae2231ee9
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_0.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── 3281ff7c3fd7441da807e780ea25565e
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_6.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── 58517f622f6c480a92f9e7c106a39057
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_3.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── 5eec59e5ed974034b25202c70462c746
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_2.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── 8880bb92c693459c8ae25d623d55f51c
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_5.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── e2d9f9fe94944261972ad263f464c19e
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_1.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   ├── ebdb493877f147c2acf4cc8660fc67bb
   │   ├── gfnff_charges
   │   ├── gfnff_topo
   │   ├── mono_8.xyz
   │   ├── output.log
   │   ├── xtbopt.log
   │   └── xtbopt.xyz
   └── htp_test.py

   9 directories, 55 files

Thus, each unique-labeled directory displayed the initial **mono_IDX.xyz** structure, and the results after the geometry optimization procedure performed 
at the **GFN-FF** level of theory. The final relaxed structure is written in the file **xtbopt.xyz**.  
   
