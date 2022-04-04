from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem import TorsionFingerprints

import MDAnalysis as mda

import numpy as np

class Mcon(object):
   """Calculates molecular conformers

   Example:
   --------


   Note
   -----
   This class requires RDKit to be installed.
   """
     
   __slots__ = 'mol','num_conf','e_max' 

   def __init__(self, mol, num_conf, e_max):
       """
       Parameters

       ----------
       mol: rdkit.Chem.rdchem.Mol
            RDKit Mol object

       num_conf: int
          The number of configurations requested to be computed.
    
       e_max: bool, optional (default False)
          If True, reports the conformer with highest energy.

       """
       self.mol=mol
       self.num_conf=num_conf
       self.e_max=e_max
     
   def conformer(self):
       """
       Calculate a molecular conformer for a given molecule. If 
       enabled the highest-energy conformer is returned. 

       Returns
       -------
       datapoint: rdkit.Chem.Mol
         RDKit Mol object

       """
       mol, num_conf=self.mol, self.num_conf 
       e_max=self.e_max

       m, energies=Memb().etkdgv3_energies(mol, num_conf)  

       # Chosing the conformer with highest energy
       max_e_index=np.argmax(np.around(energies,decimals=4))

       new_mol=Chem.Mol(mol)
       
       new_mol.RemoveAllConformers()
       new_mol.AddConformer(m.GetConformer(int(max_e_index)),
                            assignId=True)


       return new_mol
       
class Memb:
    def etkdgv3_energies(self, mol, num_conf):
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

      AllChem.EmbedMolecule(mol)
      AllChem.MMFFOptimizeMolecule(mol)
      
      m=Chem.Mol(mol)
      
      ps=molDG.ETKDGv3()
      ps.randomSeed=0xf00d
      
      cids=molDG.EmbedMultipleConfs(m,int(num_conf),ps)
      mp=AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
      AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0,
                                        mmffVariant='MMFF94s')

      energies=[]
      for cid in cids:
         ff = AllChem.MMFFGetMoleculeForceField(m, mp, confId=cid)
         energies.append(ff.CalcEnergy())

      return m, energies 
