Cloning this branch:
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

   git clone -b pysoft_docs --single-branch git@github.com:

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
    pysoft_docs tracked
  Local branch configured for 'git pull':
    pysoft_docs merges with remote pysoft_docs
  Local ref configured for 'git push':
    pysoft_docs pushes to pysoft_docs (up-to-date)

If this is **true**, then you dont need to do anything else of defining origins. Then, to push changes, you just need to do

.. code-block:: console

    git status  
    git add --all -f  (sometimes the build directory is seen in the .gitignore)
    git pull
    git push

and the changes will be done in the tutorials branch.
