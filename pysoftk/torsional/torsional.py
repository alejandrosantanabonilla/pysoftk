from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import Draw

import networkx as nx

import itertools
from itertools import chain


class Torsional(object):
    """
       A class for detecting torsional angles in a linear conjugated polymer.
          

       Examples
       --------

    """    
    def __init__(self, mol):
       """Initializes the class Torsional

       Parameters
       ----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
       """
       self.mol = mol
    
    def validate_mol(self, mol):
      """Function that validates an object as RDKit Mol object.

       Parameters
       ----------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
    
       Return
       -------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
      """
      try:
          isinstance(mol,Chem.rdchem.Mol)
      except ValueError:
          raise TypeError("{0} is an invalid molecule".format(mol))
      else:
        return mol
       
    def seek_angles(self):
      """List of torsional angles in a planar polymer.

      Return
      -------

      List: int
        List of atomic indexes for torsional angles
      """
      return self.list_dhdl_atoms()

    def __list_ring_atoms(self):
      """Function to compute the inner rings belonging to an RDKit Mol object.
       
       Return
       -------

       List : int
        List of atomic indexes for torsional angles
      """
      rdkit_mol=self.validate_mol(self.mol)
      ri=rdkit_mol.GetRingInfo()
      
      return [list(values) for idx,values in enumerate(ri.AtomRings())]

    def mol_graph(self):
      """Function to create a NetworkX graph based on an RDKit molecule.
      
      Return
      ------

      G : Networkx.Graph
        A networkx graph object updated with an RDKit Mol object.
      """
      mol=self.validate_mol(self.mol)
      G=nx.Graph()
  
      for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetFormalCharge(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   is_aromatic=atom.GetIsAromatic())

      for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bon_type=bond.GetBondType())

      return G

    def __dihedral_atoms(self):
      """Function that computes the dihedral atoms of an RDKit molecule.

      Return
      ------

      List: int
        List of atomic indexes for dihedral angles
      """
      graph = self.mol_graph()
     
      trsl_atm=[i for i in range(len(graph.nodes))
               if graph.degree[i] == 3]
 
      trs_list=[list(graph.adj[value]) for idx, value
               in enumerate(trsl_atm)]

      slc_pts=[list(set(values)-set(trsl_atm))
                   for idx, values in enumerate(trs_list)]

 
      str_end_pts= [val for sublist in slc_pts
                   for val in sublist]

      tot_comb=list(itertools.combinations(list(set(str_end_pts)),2))

      inner_tor=[tuple(set(list(graph.adj[values])) - set(trsl_atm))
                for idx, values in enumerate(trsl_atm)]

      dhdl_atms=list(set(tot_comb)-set(inner_tor))
  
      return dhdl_atms
   
    def __paths_dihedral_atoms(self):
      """Function that computes atomic indexes with 4 neighbours.

      Return
      ------

      List: int
        List of all paths produced by tuple of indexes of atoms.
      """
      graph=self.mol_graph()
      dh_atoms=self.__dihedral_atoms()
      
      dihedral_atoms=list(dh_atoms)
      
      paths=[list(nx.all_simple_paths(graph, source=tup[0],
              target=tup[1],cutoff=4)) for tup in dihedral_atoms]
  
      return paths
  
    def __comb_dihedral_atoms(self):
      """ Function to detect atoms with exact 4 neighbours.
      
      Return
      ------

      List: int
        List of atoms with 4 neigbours
      """
      paths=self.__paths_dihedral_atoms()
      
      lst=list(chain.from_iterable(paths))
      atm_comb = [values for idx,values in enumerate(lst)
                 if len(values) == 4]
    
      return atm_comb

    def __prbm_cases(self):
      """Function to detect atoms with complex environments.
      
      Return
      ------

      List: int
        List of atomic indexes for atoms with complex environment.
      """
      lst_rng_atm=self.__list_ring_atoms()
      cmb_dhr_atm=self.__comb_dihedral_atoms()
      
      l1=list([i for i in range(len(lst_rng_atm))])
      l2=list([i for i in range(len(cmb_dhr_atm))])

      idx_com=list(itertools.product(l1,l2))

      prb_cas=[tup for tup in idx_com
              if len(list(set(lst_rng_atm[tup[0]]) &
                       set(cmb_dhr_atm[tup[1]]))) >= 4]

      return prb_cas

    def list_dhdl_atoms(self):
      """Function to compute all dihedral atoms

      Return
      ------ 

      List: int
        List of atomic indexes for dihedral atoms
      """
      prb_cas=self.__prbm_cases()
      cmb_dhr_atm=self.__comb_dihedral_atoms()
      
      dlt=[cmb_dhr_atm[tup[1]] for tup in prb_cas]
 
      set_dlt=set(tuple(row) for row in dlt)
      set_ttl=set(tuple(row) for row in cmb_dhr_atm)
      
      return sorted(list(set_ttl-set_dlt))

    def mol_graph_neigh(self, mol):
      """Function to compute adjacency matrix from an RDKit Mol object.
       
       Return
       ------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
      """
      rdkit_mol=self.validate_mol(self.mol)
      Chem.Kekulize(rdkit_mol)
      
      atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
      am = Chem.GetAdjacencyMatrix(mol,useBO=True)

      for i,atom in enumerate(atoms):
        am[i,i] = atom

      G = nx.from_numpy_matrix(am)

      return G

    def show_atom_number(self, mol, label):
       """Function to set atomic numbers to an RDKit Mol object.

       Parameters
       ----------

       label : str
          Added Atomic Property to an 
          RDKit Mol Object

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
       
       Return
       ------

       mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
       """
       for atom in mol.GetAtoms():
            atom.SetProp(label, str(atom.GetIdx()))
       return mol
  
    def draw_molecule(self, mol, name_pic, atoms_idx):
        """Function to print a molecule into a PNG file.

        Parameters
        ----------

        mol : rdkit.Chem.rdchem.Mol
          RDKit Mol object
 
        name_pic : str
             Name of the created PNG object

        atoms_idx: List[tup int]
             List of tuples of atomic indexes to 
             be highligthed.
       
        Return
        -------

        None 
           rdkit.Draw.MolToFile Object
        """
        self.show_atom_number(mol, 'molAtomMapNumber')
        
        Draw.MolToFile(mol,str(name_pic)+".png",
                       includeAtomNumbers=True, highlightAtoms=atoms_idx) 
  
    def plot_trs_ang(self, name):
      """Function to plot torsional atoms from an RDKit Mol Object.
      
      Parameters
      -----------

      name : str
          Name of the RDKit Mol object 

      Return
      -------

      draw_molecule : PNG 
            A PNG RDKit Mol object
      """
      rdkit_mol=self.validate_mol(self.mol)
      lst_tor_angl=self.seek_angles()
      
      for idx, values in enumerate(lst_tor_angl):
        self.draw_molecule(rdkit_mol,str(name)+'_'+'mol_'+str(idx), values)

      print ("{} figures have been printed".format(len(lst_tor_angl)))





