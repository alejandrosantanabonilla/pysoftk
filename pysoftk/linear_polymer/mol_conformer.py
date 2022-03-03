from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem import TorsionFingerprints

import numpy as np

class Mcon(object):
   """Calculates molecular conformers

   Example:
   --------


   Note
   -----
   This class requires RDKit to be installed.
   """
     
   __slots__ = 'mol','num_conf','max_iter','e_max' 

   def __init__(self, mol, num_conf, max_iter, e_max):
       """
       Parameters

       ----------
       mol: rdkit.Chem.rdchem.Mol
            RDKit Mol object

       num_conf: int
          The number of configurations requested to be computed.

       max_iter: int
          The maximum number of iterations to be used in the FF
          relaxation.
     
       e_max: bool, optional (default False)
          If True, reports the conformer with highest energy.

       """
       self.mol=mol
       self.num_conf=num_conf
       self.max_iter=max_iter
       self.e_max=e_max
     
   def conformer(self):
       """
       Calculate molecular conformers for a given molecule. If 
       enabled the highest-energy conformer is returned. 

       Returns
       -------
       datapoint: rdkit.Chem.Mol
         RDKit Mol object

       """
       mol, num_conf=self.mol, self.num_conf, 
       max_iter, e_max=self.max_iter, self.e_max

       mol=Memb().etkdgv3(mol,num_conf)

       new_mol=Chem.Mol(mol)
       energies = AllChem.MMFFOptimizeMoleculeConfs(mol,
                                                    maxIters=int(max_iter))

       cids=[conf.GetId() for conf in mol.GetConformers()]

       energies_list=np.array([e[1] for e in energies])
       max_e_index=np.argmax(np.around(energies_list,decimals=4))
       min_e_index=np.argmin(np.around(energies_list,decimals=4))
   
       # Chosing the conformer with higher or lowest energy

       if e_max == "True":
         new_mol.AddConformer(mol.GetConformer(int(max_e_index)))

       else:
         new_mol.AddConformer(mol.GetConformer(int(min_e_index)))

       return new_mol
    
class Memb:
    def etkdgv3(self, mol, num_conf):
      """Calculate molecular configurations using 
         the RDKit-ETKDG3 method.
         
         Parameters
         ----------
         datapoint: 
        
         Returns
         -------
         datapoint: rdkit.Chem.rdistGeom.EmbedMultipleConfs
            RDKit Mol object
      """
      ps = molDG.ETKDGv3()
      ps.randomSeed = 1
      molDG.EmbedMultipleConfs(mol,int(num_conf),ps)
      return mol 
