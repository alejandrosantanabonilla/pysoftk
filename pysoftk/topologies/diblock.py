import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import copy
from .utils import *
from collections import OrderedDict

from pysoftk.linear_polymer import linear_polymer as lp
from pysoftk.linear_polymer import super_monomer  as sm
from pysoftk.linear_polymer.utils import *

class Db:
    """A class for creating a diblock copolymers from given RDKit molecules.

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

       
    def diblock_copolymer(self, len_block_A, len_block_B,
                             FF="MMFF", relax_iterations=100):

        """Function to create a diblock copolymer 


        Parameters
        -----------
           
        len_block_A: int 
            Length of the molecular block A

        len_block_B: int
            Length of the molecular block B

        FF: str
            Selected FF between MMFF or UFF

        relax_iterations: int  
            Number of iterations used for relaxing a molecular object.  

        Return
        -------

        newMol_H : rdkit.Chem.rdchem.Mol
             RDKit Mol object
        """

        ma=self.ma
        mb=self.mb
        atom=self.atom
              
        monomer = sm.Sm(ma, ma, str(atom))
        
        for i in range(len_block_A-2):
            monomer=sm.Sm(monomer.mon_to_poly(),
                                  ma, str(atom))
    
        monomer2 = sm.Sm(mb, mb, str(atom))
        
        for i in range(len_block_B-2):
             monomer2=sm.Sm(monomer2.mon_to_poly(),
                                    mb, str(atom))
    
        diblock = sm.Sm(monomer.mon_to_poly(),
                              monomer2.mon_to_poly(),
                              str(atom))

        mol=diblock.mon_to_poly()  

        newMol_H=swap_hyd(mol, relax_iterations, str(atom), FF)
        
        return newMol_H


class Pt:
    """A class for creating a Patterned polymers from a list of RDKit molecules.
    
    Examples
    ---------

    Note
    -----

    RDKit package must be installed.
    """
    
    __slots__ = ['pattern', 'mols', 'atom']

    def __init__(self, pattern, mols, atom):
       """
       Initialize this class.

       Parameters
       -----------

       pattern : str
            Variable containing the desired pattern written in a string.

       mols : list
            List containing RDKit molecular objects.

       atom : str
            Variable containing the atomic place holder for merging the molecules.
          
       """
       
       self.pattern = pattern
       self.mols = mols 
       self.atom = atom
       
    def pattern_block_poly(self, relax_iterations=100):
        """
        Function to create a polymer based on an alphabetic ordered pattern.
    

        Return
        --------

        return : mol
             A molecular pysoftk object.      

        """
        pattern= self.pattern
        mols=self.mols
        atom=self.atom

        
        names=['mol_{}+'.format(i) for i in range(1,len(mols)+1)]
        seq=pattern_mol_seq(names,pattern)
   
        numbers=['mol_{}'.format(i) for i in range(1,len(mols)+1)]
        od_mols=OrderedDict(list(zip(numbers,mols)))
 
        list_mol=[od_mols[i] for i in seq]

        outmol=Chem.CombineMols(list_mol[0],list_mol[1])

        for i in range(2,len(list_mol)):
            outmol=Chem.CombineMols(outmol,list_mol[i])

        lst_ngh=atom_neigh(outmol, str(atom))
        tpb=tuple_bonds(lst_ngh)

        proto_pol=create_pol(outmol, str(atom), tpb)

        return swap_hyd(proto_pol,relax_iterations,str(atom))


