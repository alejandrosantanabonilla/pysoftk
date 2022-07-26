from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom, rdMolTransforms
from rdkit.Chem import rdDistGeom as molDG

import numpy as np

#Disable the unnecessary RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

class Lp:
    """A class for creating a linear polymer 
       from given RdKit molecules.

       Attributes:
       -----------
       mol          The super_monomer molecule to be polarised
       atom         A place-holder atom to connect the molecule 
       n_copies     Number of copies to create the polymer
       shift        A real value to translate the super_monomer 
                    in a real-grid.

       Examples:
       ---------


       Note:
       -----
       RDKit package must be installed.
    """
    
    __slots__ = 'mol', 'atom', 'n_copies', 'shift'

    def __init__(self, mol, atom, n_copies, shift):
       """Initialize this class.
          
       Parameters
       ----------
       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
 
       atom : 'str'
          The placeholder atom to combine the molecules 
          and form a new monomer

       n_copies: 'int'
          Number of copies to be created for the provided 
          mol object. 

       shift: 'float'
          X-axis shift to translate the super_monomer
          object.
       """
       
       self.mol = mol
       self.atom = atom
       self.n_copies = n_copies
       self.shift = float(2.5) if shift is None else float(shift)

    def max_dist_mol(self):
        """Returns the maximum distance between atoms 
           from an RDKit Mol Object.

        Returns
        -------
        np.float
           maximum value from a Molecule Bound Matrix 
        """
        mol=self.mol
        bm = molDG.GetMoleculeBoundsMatrix(mol)    
        return np.amax(bm)


    def x_shift(self):
        """Re-calibrate distance between monomers
     
        Returns
        -------
        shift_final : 'float'
              Value to translates the unit repetitions 
        """
        
        shift=self.shift
        shift_final=float(self.max_dist_mol())-float(shift)
        
        return shift_final

    def _copy_mol(self):
       """Function to replicate super_monomers.

       Parameters:
       -----------
       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object

       Returns:
       --------
       fragments : `list`
          List of RDKit Mol objects
       """    
       mol=self.mol
       n_copies=self.n_copies
       
       fragments=[mol for _ in range(int(n_copies))]

       return fragments


    def _polimerisation(self, fragments):
        """Function to produce a polymer in a
           recursive manner. 

        Parameters:
        -----------
        fragments: 
            A list of RDKit objects.

        Returns:
        --------
        outmol : rdkit.Chem.rdchem.Mol
             RDKit Mol object

        """
        
        x_offset = self.x_shift()
        
        outmol=fragments[0]
        for idx, values in enumerate(fragments[1:]):
            outmol = Chem.CombineMols(outmol,values,
                                      offset=Point3D(x_offset*(idx+1),0.0, 0.0))
            
        return outmol

    def _bond_conn(self, outmol):
        """Function to peruse the bonds and connections of place-hold 
           atom within a super_monomer object.
 
        Parameters:
        ----------
        outmol : rdkit.Chem.rdchem.Mol
           RDKit Mol object

        Returns:
        --------
        Tuple : `list`
          A tuple of lists containing connections and 
          neighbours of the place-holder atom.
        
        """
        atom=self.atom
        conn=[]; neigh=[]
        
        rwmol = Chem.RWMol(outmol)

        for atoms in rwmol.GetAtoms():
            if atoms.GetSymbol() == str(atom):
               conn.append(atoms.GetIdx())
               neigh.append([(atoms.GetIdx(),nbr.GetIdx())
                             for nbr in atoms.GetNeighbors()])

        return conn[1:-1], neigh[1:-1]

    def _swap_hyd(self, mol, FF="MMFF"):
        """Function that seeks for last place-holder atoms 
           and convert them into Hydrogen atoms.

        Parameters:
        -----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit Mol object
       
        Returns:
        --------
        newMol_H : rdkit.Chem.rdchem.Mol
            RDKit Mol object

        """
       
        atom=self.atom
        
        for atoms in mol.GetAtoms():
            if atoms.GetSymbol() == str(atom):
                mol.GetAtomWithIdx(atoms.GetIdx()).SetAtomicNum(1)

        newMol = AllChem.AssignBondOrdersFromTemplate(mol, mol)
        newMol_H = Chem.AddHs(newMol, addCoords=True)

        if FF == "MMFF":
           AllChem.MMFFOptimizeMolecule(newMol_H,maxIters=5)

        else:
           AllChem.UFFOptimizeMolecule(newMol_H,maxIters=5)

        return newMol_H
   
    def linear_polymer(self, FF="MMFF", iter_ff=100):
       """Function to create a linear polymer 
          from a given super_monomer object.

       Parameters
       ------------

       FF : str
          Selected Force-Field from RDKit options, i.e, UFF or MMFF
         
       iter_ff: int
          The maximum number of iterations used for the FF.

       Returns:
       --------
       newMol_H : rdkit.Chem.rdchem.Mol
            RDKit Mol object

       """
       atom=self.atom
       mol=self.mol
       n_copies=self.n_copies
       
       fragments=self._copy_mol()
       outmol=self._polimerisation(fragments)
       new_conn,new_neigh=self._bond_conn(outmol)   

       rwmol = Chem.RWMol(outmol)
       c_bonds=[]
       
       for idx, value in enumerate(new_neigh):
           br,c=value[0]
           c_bonds.append(c)

       tuple_c_bonds=[(c_bonds[i],c_bonds[i+1])
                      for i in range(0,len(c_bonds),2)]

       for i in tuple_c_bonds:
           a,b=i
           rwmol.AddBond(a, b, Chem.BondType.SINGLE)

       for i in sorted(new_conn, key=None, reverse=True):
           rwmol.RemoveAtom(i)

       mol3=rwmol.GetMol()

       if FF == "MMFF":
           AllChem.MMFFOptimizeMolecule(mol3,maxIters=int(iter_ff))
           newMol_H=self._swap_hyd(mol3, "MMFF")

       else:
           AllChem.UFFOptimizeMolecule(mol3,maxIters=int(iter_ff))
           newMol_H=self._swap_hyd(mol3, "UFF")
       
       return newMol_H
   



