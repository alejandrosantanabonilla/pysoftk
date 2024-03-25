import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import random 
import copy

from pysoftk.linear_polymer import linear_polymer as lp
from pysoftk.linear_polymer import super_monomer  as sm

from pysoftk.tools.utils_rdkit import *
from pysoftk.tools.utils_func import *

class Rnp():
    """A class for creating a random copolymer from given RDKit molecules.
       
    Examples
    ---------


    Note
    -----

    RDKit package must be installed.
    """
    
    __slots__ = ['ma', 'mb', 'atom']

    def __init__(self, ma, mb, atom):
       """
       Initialize this class.

       Parameters
       -----------

       ma : rdkit.Chem.rdchem.Mol
            RDKit Mol object

       mb : rdkit.Chem.rdchem.Mol
            RDKit Mol object
 
       atom : str
            The placeholder atom to combine the molecules and form a new monomer
          
       """
       
       self.ma = ma
       self.mb = mb
       self.atom = atom

    def random_ab_copolymer(self, len_polymer, pA,
                            force_field="MMFF", iter_ff=100, swap_H=True):
       """ Function to build a random copolymer using an user provided 
           probability (pA) for merging the monomer ma and imposing 
           the condition pB=1-pA.


       Parameters
       -----------

       len_polymer: int
          User defined length of the polymer.

       pA: float
          User defined attaching probability of ma.

       force_field: str
           User selected FF between "MMFF" or "UFF".   

       iter_ff: int
           User defined iterations for a FF.
  
        swap_H: bool
             Indicates if the user defined atomic place holder 
             is changed to a Hydrogen atom or remain as the 
             used species.   
 
       Return
       -------    

       mol: rdkit.Chem.rdchem.Mol
          RDKit Mol object
        
       """

       ma=self.ma
       mb=self.mb
       atom=self.atom

       rand = np.random.rand()
    
       if rand<pA:
          m1 = ma
          
       else:
          m1 = mb
    
       rand = np.random.rand()
    
       if rand<pA:
          m2 = ma
       else:
          m2 = mb
        
       monomer = sm.Sm(m1, m2, str(atom))    
    
       for i in range(len_polymer-1):
          rand = np.random.rand()
    
          if rand<pA:
              m3 = ma
          else:
              m3 = mb
            
          monomer=sm.Sm(monomer.mon_to_poly(), m3, str(atom))
        
       mol=monomer.mon_to_poly()


       valid_force_fields = ("MMFF", "UFF", "MMFF94")
       if force_field not in valid_force_fields:
            raise ValueError(f"Invalid force field: {force_field}. Valid options are: {valid_force_fields}")

       # Automatically change ff if necessary:
       if force_field == "MMFF94":
           force_field = "MMFF"  # Change to default MMFF94 for MMFF

       # Validate and convert iterations and steps to integers:
       try:
           iter_ff = int(iter_ff)
           
       except ValueError:
           raise ValueError("iter_ff must be an integer.")
       
       if swap_H:
            newMol_H=swap_hyd(mol, iter_ff, str(atom), force_field)

       if not swap_H:
            newMol_H=no_swap(mol, iter_ff, force_field)

       return newMol_H
    
    def random_abc_copolymer(self, mc, len_polymer, pA,
                             pB, force_field="MMFF", iter_ff=100, swap_H=True):
        
       """ Function to build a random copolymer based on an user defined 
           probability (pA) of merging mA, pB for monomer mb, and the 
           condition pC=1-pA-pB.

       Parameters
       -----------

       mc: rdkit.Chem.rdchem.Mol
           RDKit Mol object #

       len_polymer: float
           User defined length of the polymer.

       pA: float
           User defined attaching probability of ma.

       pB: float
           User defined attaching porbability of mb.

       force_field: str
           User selected FF.

       iter_ff: int
           User defined iterations for a FF.
     
       swap_H: bool
             Indicates if the user defined atomic place holder 
             is changed to a Hydrogen atom or remain as the 
             used species.   
 
       Return
       -------
    
       mol: rdkit.Chem.rdchem.Mol
             RDKit Mol object

       """

       ma=self.ma
       mb=self.mb
       atom=self.atom
    
       rand = np.random.rand()
       pAB=pA+pB
    
       if rand<pA:
         m1 = ma

       elif rand<pAB:
         m1 = mb

       else:
         m1 = mc
    
       rand = np.random.rand()
    
       if rand<pA:
           m2 = ma

       elif rand<pAB:
           m2 = mb

       else: 
           m2 = mc
        
       monomer = sm.Sm(m1, m2, str(atom))
    
       for i in range(len_polymer-1):
          rand = np.random.rand()
    
          if rand<pA:
              m3 = ma

          elif rand<pAB:
              m3 = mb

          else:
              m3 = mc

          monomer=sm.Sm(monomer.mon_to_poly(), m3, str(atom))
       
       mol = monomer.mon_to_poly()

       valid_force_fields = ("MMFF", "UFF", "MMFF94")
       if force_field not in valid_force_fields:
            raise ValueError(f"Invalid force field: {force_field}. Valid options are: {valid_force_fields}")

       # Automatically change ff if necessary:
       if force_field == "MMFF94":
           force_field = "MMFF"  # Change to default MMFF94 for MMFF

       # Validate and convert iterations and steps to integers:
       try:
           iter_ff = int(iter_ff)
           
       except ValueError:
           raise ValueError("iter_ff must be an integer.")

       if swap_H:
           newMol_H=swap_hyd(mol, iter_ff, str(atom), force_field)

       if not swap_H:
           newMol_H=no_swap(mol, relax_iterations, force_field)

       return newMol_H
