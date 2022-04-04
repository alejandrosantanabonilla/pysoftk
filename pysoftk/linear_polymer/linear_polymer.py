from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG

import itertools
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import distances

class Lp:
    """A class for creating a linear polymer 
       from given RdKit molecules.

       Examples:
       ---------


       Note:
       -----
       RDKit and MDAnalysis packages must be installed.
    """
    
    __slots__ = 'mol', 'atom', 'n_copies'

    def __init__(self, mol, atom, n_copies):
       """Initialize this class.
          
       Parameters
       ----------
       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
 
       atom : str
          The placeholder atom to combine the molecules 
          and form a new monomer

       n_copies: int
          Number of copies to be created for the provided 
          mol object. 
       """
       self.mol = mol
       self.atom = atom
       self.n_copies = n_copies
       
    def elem(self):
        """Returns the atomic elements of 
           an RDKit Mol Object.

        Returns
        -------
        List[str]
          List of atomic symbols
        """
        mol=self.mol
        elements = [atoms.GetSymbol() for atoms in mol.GetAtoms()]
        return elements

    def atom_3d_coords(self):
        """Returns the cartesian coordinates of 
           an RDKit Mol Object.

        Returns
        -------
        List[(n,3) float]
           List of (n,3) cartesian coordinates 
        """
        mol=self.mol
        AllChem.MMFFOptimizeMolecule(mol)
        pos = [mol.GetConformer().GetAtomPosition(i)
            for i in range(mol.GetNumAtoms())]
        return pos

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

    def vdw_val(self):
        """Returns the vdW radius of a given atom
           from an RDKit Mol Object.
        
        Returns
        -------
        dict : dict
          Dictionary of {str: float}
          with key:str as an atomic element and item:float 
          as vdw radii. 
        """
        from rdkit.Chem import PeriodicTable as pt
        atom=self.atom
        atomic_num=Chem.MolFromSmiles(str(atom)).GetAtomWithIdx(0).GetAtomicNum()
        vdw_radii=pt.GetRvdw(Chem.GetPeriodicTable(),int(atomic_num))
        return {str(atom):float(vdw_radii)}
    
    def grid_translation(self):
        """Returns a tuple of cartesian coordinates and elements
           Grid translated along x-axis in real space from an 
           RDKit Mol Object.

        Returns
        -------
        List[(n,3) float]
           List of (n,3) cartesian coordinates

        List[n str]
           List of n atomic elements 
        """
        n_copies=self.n_copies
        pos = self.atom_3d_coords()
        symbols = self.elem()

        n_atoms = n_copies * len(symbols)
        elements=symbols*n_copies
        molecule = np.array([[pos[i].x,pos[i].y,pos[i].z]
                              for i in range(len(pos))])

        max_dist=self.max_dist_mol()
        grid_size = n_copies*(max_dist+0.75) 
        spacing = max_dist+3.25

        # translating coordinates around a grid
        coordinates = []

        for i in range(n_copies):
          x = spacing * (i % grid_size)
          xyz = np.array([x, 0, 0])
          coordinates.extend(molecule + xyz.T)

        return coordinates, elements

    def mda_creation(self):
       """Create a MDAnalysis molecule object.

       Returns
       -------
       mda.pol:: MDAnalysis.mda
            A MDAnalysis object containing the information of 
            a molecule.
       """
       res_1, res_2=self.grid_translation()
       pos_end=len(res_1)
       mda_pol= mda.Universe.empty(pos_end,
                               trajectory=True)
 
       mda_pol.atoms.positions=np.array(res_1)
       mda_pol.add_TopologyAttr('names',res_2)
       mda_pol.add_TopologyAttr('elements',res_2)
       mda_pol.add_TopologyAttr('types',res_2)
       return mda_pol
       
    def mda_polymer(self):
       """
       Returns
       -------
       mda.pol:: MDAnalysis.mda
            The constructed molecule as a MDAnalysis object.

       brt_ext:: List[tup int]
            List of tuples with indexes of atom placeholders.
       """
       atom=self.atom
       vdw_list=self.vdw_val()
       mda_pol=self.mda_creation()
       
       query=" ".join(["name",str(atom)])
       br_dist=mda_pol.select_atoms(str(query))
       
       tmp_br=distances.self_distance_array(br_dist.positions)
       dist_br=np.around(tmp_br, decimals=2)

       br_atoms=[i.index for i in br_dist.atoms]
       br_comb=list(itertools.combinations(br_atoms,2))

       result=np.transpose(np.where(dist_br == np.amin(dist_br)))
       result_idx= [i[0] for i in result]

       bond_br=[br_comb[i] for i in result_idx]

       cc_bonds=[(values[0]-1,values[1]-1)
                     for idx, values in enumerate(bond_br)]
   
       bonds_core=mda.topology.guessers.guess_bonds(mda_pol.atoms,
                                                    mda_pol.atoms.positions,
                                                    vdwradii=vdw_list)

       max_separation=np.transpose(np.where(dist_br == np.amax(dist_br)))
       result_max_br = [i[0] for i in max_separation]
       bond_max_br=[br_comb[i] for i in result_max_br]
       br_ext = [tup for tup in bond_max_br[0]]

       bonds_total=cc_bonds + list(bonds_core)
       mda_pol.add_TopologyAttr('bonds',bonds_total)
       
       return mda_pol, br_ext
   
    def detect_mda_polymer(self): 
       """ Create a linear polymer 
      
       Returns
       -------
       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
       """
       atom=self.atom
       result,br_ext=self.mda_polymer()
       
       mol1=result.atoms.convert_to("RDKIT",NoImplicit=False)

       all_br=[atoms.GetIdx() for atoms in mol1.GetAtoms()
                          if atoms.GetSymbol() == str(atom)]

       br_rmv=sorted(list(set(all_br)-set(br_ext)),reverse=True)

       emol = Chem.EditableMol(mol1)
       for atom in br_rmv :
           emol.RemoveAtom(atom)

       mol2=emol.GetMol()
       return mol2

    def linear_polymer(self,ff_itrs=100): 
        """Construct a linear polymer.  

        Parameters
        ----------
        fft_itrs : int
          RDKit Force Field iteration cycles.     

        Returns
        -------
        mol : rdkit.Chem.rdchem.Mol
           RDKit Mol object
        """
        atom=self.atom
        mol_con=self.detect_mda_polymer()
        
        mol=Chem.RWMol(mol_con)

        for atoms in mol.GetAtoms():
            if atoms.GetSymbol() == str(atom):
               mol.GetAtomWithIdx(atoms.GetIdx()).SetAtomicNum(1)

        mol2=mol.GetMol()

        AllChem.UFFOptimizeMolecule(mol2,maxIters=int(ff_itrs))

        return mol2




  



   


















