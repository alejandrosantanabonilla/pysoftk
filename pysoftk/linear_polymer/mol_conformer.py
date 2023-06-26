from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem import TorsionFingerprints

from pysoftk.tools.utils_rdkit import *
import numpy as np

class Mcon(object):
   """Calculates molecular conformers

   Example
   --------


   Note
   -----

   This class requires RDKit to be installed.
   """
     
   __slots__ = ['mol','num_conf'] 

   def __init__(self, mol, num_conf):
       """
       Parameters
       ----------

       mol: rdkit.Chem.rdchem.Mol
            RDKit Mol object

       num_conf: int
          The number of configurations requested to be computed.
    
       """
       self.mol=mol
       self.num_conf=num_conf
 
   def conformer(self, output_name):
       """
       Calculate a molecular conformer for a given molecule. If 
       enabled the highest-energy conformer is returned. 

       Return
       -------

       list of molecules: rdkit.Chem.Mol
         SDWriter Mol object

       """
       mol=self.mol
       num_conf=self.num_conf 

       mol, cids, energies=etkdgv3_energies(mol, num_conf)  

       full_name=[output_name,".sdf"]
       w=AllChem.SDWriter("".join(full_name)) 

       for cid in cids:
          w.write(mol, confId=cid)   
       w.close()
       
       print ("Number of conformers created: {} ".format(len(cids)))
