from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG

#Disable the unnecessary RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

class Sm(object):
    """ Class to create a new combined molecule from two provided RDKit Mol objects.

    Example
    -------


    Note
    ----
    This class requires RDKit to be installed.
    """
    
    __slots__ = ['mol_1', 'mol_2' , 'atom']

    def __init__(self, mol_1, mol_2, atom):
       """Initializes the class Sm.
          
       Parameters
       ----------

       mol_1 : rdkit.Chem.rdchem.Mol
          First molecule for the new monomer

       mol_2 : rdkit.Chem.rdchem.Mol
          Second molecule for the new monomer
       
       atom : `str`
          The placeholder atom to combine the molecules and form a new monomer
       """
       self.mol_1 = mol_1
       self.mol_2 = mol_2
       self.atom = atom
       
    def constructor(self):
        """Function to combine two molecules using an atom placeholder.
        
        Returns
        -------

        mol_4 : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol Object.
        """
        mol_1, mol_2, atom = self.mol_1, self.mol_2, self.atom

        patt=['[*:1]','[',str(atom),']','.','[*:2]',
        '[',str(atom),']','>>','[*:1]','[*:2]']
        react="".join(patt)
        rxn = AllChem.ReactionFromSmarts(str(react))

        m3=[]
        results = rxn.RunReactants([mol_1, mol_2])

        for products in results:
            for mol in products:
                m3.append(Chem.MolToSmiles(mol))

        mol4=Chem.MolFromSmiles(m3[0])

        return mol4

    def _bond_order(self, mol):
       """Function to create bond template for Hydrogen atoms addition.

       Parameters
       ----------

       mol : rdkit.Chem.rdchem.Mol
          First molecule for the new monomer
       
       Returns
       -------

       mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol Object.
       
       """
       
       newMol = AllChem.AssignBondOrdersFromTemplate(mol, mol)
       newMol_H = Chem.AddHs(newMol)
       
       return newMol_H
       
    def mon_to_poly(self):
       """ Function to create a RDKit Mol object to prepare for x-axis translation.

       Returns
       -------

       mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol Object.
       """
       mol=self.constructor()
       newMol_H=self._bond_order(mol)
       AllChem.EmbedMolecule(newMol_H)
       #Emb().etkdgv3(newMol_H)
       return newMol_H

    def monomer(self):
        """ Function to produce a single polymer unit as an RDKit Mol Object.
        
       Returns
       -------

          mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol Object.
        """
        atom=self.atom
        mol=Chem.RWMol(self.mon_to_poly())
        
        for atoms in mol.GetAtoms():
            if atoms.GetSymbol() == str(atom):
               mol.GetAtomWithIdx(atoms.GetIdx()).SetAtomicNum(1)

        mol2=mol.GetMol()
        newMol_H=self._bond_order(mol2)
        AllChem.EmbedMolecule(newMol_H)
        #Emb().etkdgv3(newMol_H)

        return newMol_H

   
class Emb:
    """ Class created to use the ETKDGV3 method provided in RDKit for 3D Geometry embedding.  
    """
    def etkdgv3(self, mol):
      """Function to perform Geometry embedding.

      Parameters
      ----------

      mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol Object.

      Returns
      -------

      mol : rdkit.Chem.rdchem.Mol
          The resulting embedded molecule as RDKit Mol Object.
      """
      ps=molDG.ETKDGv3()
      ps.useExpTorsionAnglePrefs = True
      ps.useBasicKnowledge = True
      ps.enforceChirality = True
      ps.randomSeed = 1
      
      molDG.EmbedMolecule(mol,ps)

      return mol 

