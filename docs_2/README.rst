
PySoftK Documentation
========================

In this folder, a copy of the source docs is stored as part of the PySoftK development workflow. This can be used to update the future 
version of the documentation as follows:

Compiling the HTML version:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To compile and get locally the Documents, these are the commands:

.. code-block:: console

   [~]$ make clean
   [~]$ make html

The **html** results can be displayed via an internet browser (chrome, firefox, etc), following these lines of code

.. code-block:: console

   [~]$ cd build/html 
   [~]$ firefox intro.html


Cloning this branch:
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

   git clone -b pysoftk_docs --single-branch https://github.com/alejandrosantanabonilla/pysoftk.git

Pushing into this branch:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

make changes, if you have used the -b tutorials then 

.. code-block:: console

   git remote show origin 

must show the following result

.. code-block:: console

   * remote origin
  Fetch URL: https://github.com/alejandrosantanabonilla/pysoftk.git
  Push  URL: https://github.com/alejandrosantanabonilla/pysoftk.git
  HEAD branch: master
  Remote branch:
    tutorials tracked
  Local branch configured for 'git pull':
    tutorials merges with remote tutorials
  Local ref configured for 'git push':
    tutorials pushes to tutorials (up-to-date)

If this is **true**, then you dont need to do anything else of defining origins. Then, to push changes, you just need to do

.. code-block:: console

    git status  
    git add --all -f  (sometimes the build directory is seen in the .gitignore)
    git pull
    git push

and the changes will be done in the tutorials branch.


