import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import copy
from collections import OrderedDict

from pysoftk.linear_polymer import linear_polymer as lp
from pysoftk.linear_polymer import super_monomer  as sm

from pysoftk.tools.utils_rdkit import *
from pysoftk.tools.utils_func import *

from pysoftk.topologies.ring import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

class Db:
    """A class for creating a diblock copolymers from 
       given RDKit molecules.

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
            The placeholder atom to combine the molecules 
            and form a new monomer
       
          
       """
       
       self.ma = ma
       self.mb = mb
       self.atom = atom 

       
    def diblock_copolymer(self, len_block_A, len_block_B,
                             force_field="MMFF", relax_iterations=100, rot_steps=1):

        """Function to create a diblock copolymer 


        Parameters
        -----------
           
        len_block_A: int 
            Length of the molecular block A

        len_block_B: int
            Length of the molecular block B

        force_field: str
            Selected force field between MMFF or UFF

        relax_iterations: int  
            Number of iterations used for relaxing a 
            molecular object.

        rot_steps: int
             Number of rotational steps for finding conformations.

        Return
        -------

        diblock : rdkit.Chem.rdchem.Mol
             RDKit Mol object
        """

        ma=self.ma
        mb=self.mb
        atom=self.atom

        string_1='A'*len_block_A
        monomer=Pt(string_1, [ma], str(atom)).pattern_block_poly(relax_iterations, force_field=str(force_field), swap_H=False)
        

        string_2='B'*len_block_B
        monomer2=Pt(string_2, [mb], str(atom)).pattern_block_poly(relax_iterations, force_field=str(force_field), swap_H=False)

        mon_temp=monomer.write("mol")
        mon2_temp=monomer2.write("mol")

        mon=Chem.MolFromMolBlock(mon_temp)
        mon2=Chem.MolFromMolBlock(mon2_temp)
        
        diblock=sm.Sm(mon, mon2, str(atom)).monomer()

        last_mol=check_proto(diblock, force_field, relax_iterations, rot_steps)
            
        return last_mol


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
       
    def pattern_block_poly(self, relax_iterations=100, force_field="MMFF", swap_H=True, rot_steps=1, more_iter=10):
        """
        Function to create a polymer based on an alphabetic ordered pattern.
    
        Parameters
        -----------

        relax_iterations: int  
            Number of iterations used for relaxing a molecular object.

        force_field: str
            Selected FF between MMFF or UFF

        swap_H: bool
             Indicates if the user defined atomic place holder is changed to a 
             Hydrogen atom or remain as the used species. 

        rot_steps: int
             Number of rotational steps for finding conformations.

        more_iter: int
             Number of extra iterations for further optimising the molecular 
             object.

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
        
        if swap_H:
            newMol_H=swap_hyd(proto_pol, relax_iterations, str(atom), force_field)

        if not swap_H:
            newMol_H=no_swap(proto_pol, relax_iterations, force_field)


        last_mol=check_proto(newMol_H, force_field, relax_iterations*more_iter, rot_steps)
            
        return last_mol
