PySoftK
=============

PysoftK is a set of Python tools and programs for modelling and simulating polymers with different topologies. The program is still under active 
development and contributions are welcome. To install pysoftk, we encourage to do it inside a virtual environment, which can be achieved in the following 
way:

1.) Create a directory named as you want and access it (in this case called work_pol):

.. code-block:: bash
 
  mkdir work_pol
  cd work_pol

2.) Create a virtual environment named pol (this name can be change, of course)

.. code-block:: bash

   python -m venv pol
   source pol/bin/activate

3.) Get pysfotk from the GitHub repository

.. code-block:: bash

  git clone https://github.com/alejandrosantanabonilla/pysoftk.git


4.) Install pysoftk in this folder using the virtual environment

.. code-block:: bash

  pip install -e .

5.) Testing pysoftk

You need to go to the folder test and then type:

.. code-block:: bash

  pytest




