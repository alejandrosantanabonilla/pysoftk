.. include:: include/links.rst

Usage
=====

.. _installation:

Installation
------------

To install PySoftK_, we encourage to do it inside a virtual environment, which can be achieved in the following way.

1. Created a directory, for instance work_pol
   
.. code-block:: bash
   
   $ mkdir work_pol

2. Get inside the created directory (work_pol)

.. code-block:: bash

   $ cd work_pol

3. Create a virtual environment name polymer (as an example)

.. code-block:: bash

   $ python -m venv polymer
   $ source polymer/bin/activate

4. Download PySoftK_ from GitHub_:

.. code-block:: console

   (.venv) $ git clone https://github.com/alejandrosantanabonilla/pysoftk.git
   
5. Install PySoftK_ using pip:

.. code-block:: console

   (.venv) $ pip install -e .
   

6. Testing PySoftK_ can be done using pytest within the test folder:

.. code-block:: console

   $ cd test
   (.venv) $ pytest
