Welcome to PySoftK's documentation!
=======================================

Python Soft-matter King's College London (PySoftK_) is a set of tools and Python_ modules for setting up, manipulating, and running atomistic simulations of
polymers. The code is freely available under the GNU_ General Public License
v3.0. 

.. _Python: https://www.python.org/
.. _PySoftK: https://github.com/alejandrosantanabonilla/pysoftk
.. _GNU: https://www.gnu.org/licenses/gpl-3.0.en.html

PySoftK_ provides interfaces to different codes through :mod:`pysoftk.linear\_polymer.calculators` which are used together with :mod:`pysoftk.linear\_polymer.super\_monomer` for creating and running atomistic molecular dynamics simulations for a broad range of polymeric materials.  

By simply inputting the SMILES code(s) for the desired monomer(s), the code allows the user to build polymers with different compositions (homopolymers, block heterpolymers, sequenced heteropolymers and random heteropolymers) and topologies (linear and ring).  In the end the code generates an equilibrated structure of the polymer that can be used as an input for building the desired initial configuration for a molecular dynamcis simulation.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Contents:
---------

.. toctree::

   about
   usage
   tutorials/tutorials
   faq
   modules
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

  


  

