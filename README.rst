PySoftK
=============

PysoftK is a set of Python tools and programs for modelling and simulating polymers with different topologies. The program is still under active 
development and contributions are welcome. To install pysoftk, we encourage to do it inside a virtual environment, which can be achieved in the following 
way:

1.) Create a directory named as you want, in this case, we call it, work_pol

mkdir work_pol

then type

cd work_pol

2.) Create a virtual environment named pol (this name can be change, of course)

.. code_block:: console

   python -m venv pol
   source pol/bin/activate

3.) Get pysfotk from the GitHub repository

.. code_block:: console

  git clone https://github.com/alejandrosantanabonilla/pysoftk.git


4.) Install pysoftk in this folder using the virtual environment

.. code_bloack:: console

  pip install -e .

5.) Testing pysoftk

You need to go to the folder test and then follow the instructions there.
