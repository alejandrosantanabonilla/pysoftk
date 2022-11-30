.. include:: /include/links.rst

.. _folders:

=================================================================
Automatic Folder Creation
=================================================================

The module :mod:`pysoftk.folder\_manager.folder\_creator` is a set of Python_ tools that provides a range of functionalities to carry out high-throughput calculations (HTPc). 

In the first instance, the creation of a user-defined number of folders with unique labels can be achieved in the following lines of code: 

.. literalinclude:: scripts/test_folder.py
   :lines: 1-5

The function receives as the first argument the number of folders to be created, whilst the second argument declares the number of cores used to perform this operation. We can verify the creation of these folders by typing the following command:

.. code-block:: bash
   
   $ tree
   .
   ├── 37c49461ba7b4206ba5cc5c101485faa
   ├── f819d762c1234f89a1bd762022a5007c
   └── test_fol_m.py


where it can be seen that the two folders (with unique alphanumeric labels) have been succesfully created in a folder where the **test_fol_m.py** is located.
   
Similarly, PySoftK_ is able to seek for files with a user-provided extension (such as .smi, .pdb, .mol, .xyz, etc). 

.. literalinclude:: scripts/test_folder.py
   :lines: 6-14

The variable **a** stores the *absolute paths* where the files **mol_1.smi** and **mol_2.smi** are stored. Likewise, the length of this list indicates the number of files with that given extension present in the current working directory.

Complementing these previous functions, PySoftK_ enables the automatic moving of the files into unique labeled folders as indicated in the following snipet:


.. literalinclude:: scripts/test_folder.py
   :lines: 15-22

To verify that the files have correctly moved, we can use the **tree** bash command as indicated above:

.. code-block:: bash

   $ tree
   
   ├── 8fa730a5245c4656ac305a7ecc9bba07
   │   └── mol_2.smi
   ├── b38be2d4de134b49915c7d8c694959f3
   │   └── mol_1.smi
   └── test_folder.py

   2 directories, 3 files

As shown above, the files **mol_1.smi** and **mol_2.smi** have been correcty allocated into unique folders. This collection of PySoftK_ tools has been 
developed to be used independent of the Polymer module and intended to be used as stand-alone or in conjunction with the HTPc tools.
