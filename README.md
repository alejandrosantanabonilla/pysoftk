<!-- Pytest Coverage Comment:Begin -->
\n
<!-- Pytest Coverage Comment:End -->

PySoftK
=============

**PySoftK** is a set of Python tools and programs for modelling and simulating polymers with different topologies. The program is still under active 
development and contributions are welcome. A complete introduction into the program can be found in this link Documentation_. To quickly install **PySoftk**, we encourage to do it inside a virtual environment, which can be achieved in the 
following way:

1.) Create a directory named as you want and access it (in this case called work_pol):

.. code-block:: bash
 
  [~] mkdir work_pol
  [~] cd work_pol

2.) Create a virtual environment named pol (this name can be changed, of course) and activate the environment:

.. code-block:: bash

   [~] python -m venv pol
   [~] source pol/bin/activate

3.) Get PySoftK from the GitHub repository:

.. code-block:: bash

  git clone https://github.com/alejandrosantanabonilla/pysoftk.git


4.) Install PySoftK in this folder using the virtual environment

.. code-block:: bash

   [~] cd pysoftk
   [~] pip install .

**NOTE:** To use the **calculators** option, the code xtb_ needs to be installed and linked to the executable using the command:

.. code-block:: bash
 
   export PATH="$PATH:$HOME/PATH_TO_FOLDER/xtb-X.X.X/bin"
   export XTBHOME="PATH_TO_XTB_FOLDER/xtb-X.X.X"
   ulimit -s unlimited

Replacing **PATH_TO_FOLDER** by the actual path where xtb is stored in your computer.

  
5.) For testing PySoftK, you need to go to the folder **test** and then type:

.. code-block:: bash

  pytest


.. _Documentation: https://alejandrosantanabonilla.github.io/pysoftk/
.. _xtb: https://github.com/grimme-lab/xtb
