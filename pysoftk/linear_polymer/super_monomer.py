from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG

class Sm(object):
    """ Class to create a new combined molecule
        from two provided RDKit Mol objects.

    Example
    -------


    Note
    ----
    This class requires RDKit to be installed.
    """
    
    __slots__ = 'mol_1', 'mol_2' , 'atom'

    def __init__(self, mol_1, mol_2, atom):
       """Initializes the class Sm.
          
       Parameters
       ----------
       mol_1 : rdkit.Chem.rdchem.Mol
          First molecule for the new monomer

       mol_2 : rdkit.Chem.rdchem.Mol
          Second molecule for the new monomer
       
       atom : `str`
          The placeholder atom to combine the molecules 
          and form a new monomer
       """
       self.mol_1 = mol_1
       self.mol_2 = mol_2
       self.atom = atom
       
    def constructor(self):
        """Function to combine two molecules using an atom
           placeholder.
        
        Returns
        -------
        mol_4 : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol
          Object.
        """
        mol_1, mol_2, atom = self.mol_1, self.mol_2, self.atom
        combined = Chem.CombineMols(mol_1, mol_2)
        rwmol = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(combined)))
        conn=[]
        neigh=[]
        for atoms in rwmol.GetAtoms():
           if atoms.GetSymbol() == str(atom):
               conn.append(atoms.GetIdx())
               neigh.append([(atoms.GetIdx(),nbr.GetIdx())
                       for nbr in atoms.GetNeighbors()])
    
        br_1,c_1=neigh[0][0]
        br_3,c_3=neigh[2][0]
    
        rwmol.AddBond(c_1, c_3, Chem.BondType.SINGLE)
    
        for i in range(0,len(conn),2):
            rwmol.GetAtomWithIdx(conn[i]).SetAtomicNum(0)
    
        mol3=Chem.DeleteSubstructs(rwmol, Chem.MolFromSmarts('[#0]'))
        mol4=Chem.AddHs(mol3)
        
        return mol4
    
    def mon_to_poly(self):
       """ Function to create a RDKit Mol object to
           prepare for grid translation.

       Returns
       -------
       mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol
          Object.
       """
       #mol_1, mol_2, atom = self.mol_1, self.mol_2, self.atom
       mol=self.constructor()#mol_1, mol_2, atom)
       Emb().etkdgv3(mol)
       return mol

    def monomer(self):
        """ Function to produce a single polymer unit
            as an RDKit Mol Object.
        
       Returns
       -------
          mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol
          Object.
        """
        atom=self.atom
        mol=Chem.RWMol(self.mon_to_poly())
        
        for atoms in mol.GetAtoms():
            if atoms.GetSymbol() == str(atom):
               mol.GetAtomWithIdx(atoms.GetIdx()).SetAtomicNum(1)

        Chem.SanitizeMol(mol)
        Emb().etkdgv3(mol)

        return mol

   
class Emb:
    """ Class created to use the ETKDGV3 
        method provided in RDKit for 3D
        Geometry embedding.  
    """
    def etkdgv3(self, mol):
      """Function to perform Geometry embedding.

      Parameters
      ----------
      mol : rdkit.Chem.rdchem.Mol
          The resulting combined molecule as RDKit Mol
          Object.

      Returns
      -------
      mol : rdkit.Chem.rdchem.Mol
          The resulting embedded molecule as RDKit Mol
          Object.
      """
      ps=molDG.ETKDGv3()
      ps.useExpTorsionAnglePrefs = True
      ps.useBasicKnowledge = True
      ps.enforceChirality = True
      ps.randomSeed = 1
      
      molDG.EmbedMolecule(mol,ps)

      return mol 















