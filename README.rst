Cloning this branch:
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

   git clone -b  --single-branch git@github.com:

Pushing into this branch:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

make changes, if you have used the -b tutorials then 

.. code-block:: console

   git remote show origin 

must show the following result

.. code-block:: console

   * remote origin
  Fetch URL: git@github.com:
  Push  URL: git@github.com:
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
