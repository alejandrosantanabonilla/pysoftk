from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import Draw

import networkx as nx

import itertools
from itertools import chain


class Torsional(object):
    """
    A class dedicated to **detecting and visualizing torsional angles** within
    molecular structures, particularly useful for understanding the conformational
    flexibility of linear conjugated polymers.

    This class provides tools to:
    1. Validate RDKit molecular objects.
    2. Convert RDKit molecules into NetworkX graphs for topological analysis.
    3. Identify and list atoms involved in rotatable bonds, which contribute to
       torsional angles.
    4. Handle complex ring systems to accurately identify relevant torsional angles.
    5. Visualize the molecule with highlighted torsional angles by generating PNG images.

    Examples
    --------
    >>> from rdkit import Chem
    >>> from pysoftk.torsional.torsional import Torsional
    >>>
    >>> # Create a simple molecule (e.g., a short polymer segment)
    >>> mol = Chem.MolFromSmiles("c1ccccc1-c2ccccc2-c3ccccc3") # Biphenyl trimer
    >>> mol = Chem.AddHs(mol) # Add hydrogens for correct atom degrees
    >>>
    >>> # Initialize the Torsional class
    >>> torsional_analyzer = Torsional(mol=mol)
    >>>
    >>> # Get the list of atoms involved in torsional angles
    >>> torsional_atoms = torsional_analyzer.seek_angles()
    >>> print(f"Identified torsional angles (atom indices): {torsional_atoms}")
    >>>
    >>> # Plot and save images highlighting these torsional angles
    >>> torsional_analyzer.plot_trs_ang(name="biphenyl_trimer_torsions")
    3 figures have been printed (example output)

    Note
    -----
    This class **requires the RDKit and NetworkX packages to be installed**
    in your Python environment. It is particularly designed for linear-like
    polymer structures where the concept of "torsional angles" between repeating
    units or significant backbone bonds is relevant.
    """
    
    def __init__(self, mol):
        """
        Initializes the `Torsional` class with the RDKit molecule to be analyzed.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the molecule for which
            torsional angles need to be detected. It is recommended to
            add hydrogens to the molecule (`Chem.AddHs(mol)`) before
            initialization for accurate bond perception and degree calculation.
        """
        
        self.mol = mol

    def validate_mol(self, mol):
        """
        Validates whether the input object is a valid RDKit Mol object.

        This internal helper function ensures that subsequent operations
        are performed on a proper RDKit molecular structure, preventing
        potential errors.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The object to be validated as an RDKit Mol object.

        Return
        -------
        mol : rdkit.Chem.rdchem.Mol
            The **validated RDKit Mol object** if the input is valid.

        Raises
        ------
        TypeError
            If the input `mol` is not an instance of `rdkit.Chem.rdchem.Mol`.
        """
        
        try:
            isinstance(mol,Chem.rdchem.Mol)
        except ValueError:
            raise TypeError("{0} is an invalid molecule".format(mol))
        else:
            return mol

    def seek_angles(self):
        """
        Identifies and returns a list of atom indices that define the
        torsional angles within the molecule, specifically focusing on
        rotatable bonds in a planar or linear polymer context.

        This is the primary public method to obtain the identified torsional angles.
        It orchestrates several internal methods to first identify potential
        dihedral atoms, then filter out those involved in problematic (e.g., ring)
        environments, to provide a clean list of relevant torsional angles.

        Return
        -------
        List : List[int]
            A **list of lists**, where each inner list contains four integer
            atom indices `[i, j, k, l]` that define a torsional angle. The
            torsional angle is the angle between the planes defined by
            atoms `(i, j, k)` and `(j, k, l)`.
        """
        
        return self.list_dhdl_atoms()

    def __list_ring_atoms(self):
        """
        Computes and returns a list of atom indices belonging to the rings
        within the RDKit molecule.

        This internal method helps identify cyclic parts of the molecule,
        which are crucial for distinguishing between rotatable bonds in linear
        chains and fixed orientations within rings.

        Return
        -------
        List : List[List[int]]
            A **list of lists**, where each inner list contains the integer
            atom indices of atoms forming a single ring in the molecule.
        """
        
        rdkit_mol = self.validate_mol(self.mol)
        ri = rdkit_mol.GetRingInfo()

        return [list(values) for idx, values in enumerate(ri.AtomRings())]

    def mol_graph(self):
        """
        Converts the RDKit molecule into a NetworkX graph representation.

        This method facilitates graph-based algorithms for analyzing molecular
        topology, such as finding paths and atom degrees, which are essential
        for identifying dihedral angles. Atom and bond properties from RDKit
        are stored as node and edge attributes in the NetworkX graph.

        Return
        ------
        G : networkx.Graph
            A **NetworkX graph object** where nodes represent atoms (with
            attributes like atomic number, formal charge, etc.) and edges
            represent bonds (with attributes like bond type).
        """
        
        mol = self.validate_mol(self.mol)
        G = nx.Graph()

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
        """
        Computes the candidate sets of four atom indices that could define
        dihedral angles in the molecule.

        This internal method identifies potential dihedral angles by looking
        for atoms with specific degrees and their neighbors within the molecular graph.

        Return
        ------
        List : List[Tuple[int, int]]
            A **list of tuples**, where each tuple contains two integer atom
            indices `(a, b)` that are potential "endpoints" of a dihedral
            angle path.
        """
        
        graph = self.mol_graph()

        trsl_atm = [i for i in range(len(graph.nodes))
                    if graph.degree[i] == 3]

        trs_list = [list(graph.adj[value]) for idx, value
                    in enumerate(trsl_atm)]

        slc_pts = [list(set(values) - set(trsl_atm))
                   for idx, values in enumerate(trs_list)]

        str_end_pts = [val for sublist in slc_pts
                       for val in sublist]

        tot_comb = list(itertools.combinations(list(set(str_end_pts)), 2))

        inner_tor = [tuple(set(list(graph.adj[values])) - set(trsl_atm))
                     for idx, values in enumerate(trsl_atm)]

        dhdl_atms = list(set(tot_comb) - set(inner_tor))

        return dhdl_atms

    def __paths_dihedral_atoms(self):
        """
        Computes all simple paths of length 4 (i.e., involving 4 atoms for a dihedral angle)
        between the identified dihedral atom pairs.

        This internal method expands on the `__dihedral_atoms` results by
        finding the actual 4-atom sequences that define dihedral angles, ensuring
        they form a continuous chain.

        Return
        ------
        List : List[List[int]]
            A **list of lists**, where each inner list represents a path of
            four integer atom indices `[atom1, atom2, atom3, atom4]` that
            constitute a potential torsional angle.
        """
        
        graph = self.mol_graph()
        dh_atoms = self.__dihedral_atoms()

        dihedral_atoms = list(dh_atoms)

        paths = [list(nx.all_simple_paths(graph, source=tup[0],
                                          target=tup[1], cutoff=4)) for tup in dihedral_atoms]

        return paths

    def __comb_dihedral_atoms(self):
        """
        Filters the paths to identify only those that are exactly 4 atoms long,
        which directly correspond to the atom sequences defining dihedral angles.

        Return
        ------
        List : List[List[int]]
            A **list of lists**, where each inner list contains exactly four
            integer atom indices that define a specific dihedral angle.
        """
        
        paths = self.__paths_dihedral_atoms()

        lst = list(chain.from_iterable(paths))
        atm_comb = [values for idx, values in enumerate(lst)
                    if len(values) == 4]

        return atm_comb

    def __prbm_cases(self):
        """
        Identifies "problematic" cases where potential dihedral angles are
        part of a ring system, meaning they do not represent a freely rotatable
        bond in the typical sense.

        This internal method is critical for filtering out non-rotatable
        "dihedral" angles that are constrained by ring closures, ensuring
        only true torsional angles are reported.

        Return
        ------
        List : List[Tuple[int, int]]
            A **list of tuples**, where each tuple `(idx_ring, idx_dihedral)`
            indicates that the dihedral angle identified by `cmb_dhr_atm[idx_dihedral]`
            has at least four of its atoms (or all its atoms) also present in
            the ring system `lst_rng_atm[idx_ring]`. These are cases to be
            excluded from the final list of rotatable torsional angles.
        """
        
        lst_rng_atm = self.__list_ring_atoms()
        cmb_dhr_atm = self.__comb_dihedral_atoms()

        l1 = list([i for i in range(len(lst_rng_atm))])
        l2 = list([i for i in range(len(cmb_dhr_atm))])

        idx_com = list(itertools.product(l1, l2))

        prb_cas = [tup for tup in idx_com
                   if len(list(set(lst_rng_atm[tup[0]]) &
                               set(cmb_dhr_atm[tup[1]]))) >= 4]

        return prb_cas

    def list_dhdl_atoms(self):
        """
        Computes and returns the final list of valid torsional angles
        (defined by four atom indices) by excluding those that are part
        of ring systems.

        This method combines the results from `__prbm_cases` and
        `__comb_dihedral_atoms` to provide the definitive list of
        rotatable torsional angles in the molecule.

        Return
        ------
        List : List[List[int]]
            A **sorted list of lists**, where each inner list contains four
            integer atom indices representing a distinct rotatable torsional angle.
            These angles are typically associated with bonds that can undergo
            conformational changes.
        """
        
        prb_cas = self.__prbm_cases()
        cmb_dhr_atm = self.__comb_dihedral_atoms()

        dlt = [cmb_dhr_atm[tup[1]] for tup in prb_cas]

        set_dlt = set(tuple(row) for row in dlt)
        set_ttl = set(tuple(row) for row in cmb_dhr_atm)

        return sorted(list(set_ttl - set_dlt))

    def mol_graph_neigh(self, mol):
        """
        Computes the adjacency matrix (with atomic numbers on the diagonal)
        from an RDKit Mol object and converts it into a NetworkX graph.

        This function creates a graph representation where nodes are atoms
        and edges represent bonds, useful for various graph-based analyses,
        including those involving connectivity and neighborhoods.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The RDKit Mol object for which the adjacency matrix and graph
            are to be computed.

        Return
        ------
        G : networkx.Graph
            A **NetworkX graph object** representing the molecular connectivity.
            Node attributes include atomic numbers (on the diagonal of the adjacency
            matrix used to construct the graph).
        """
        
        rdkit_mol = self.validate_mol(self.mol)
        Chem.Kekulize(rdkit_mol)

        atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        am = Chem.GetAdjacencyMatrix(mol, useBO=True)

        for i, atom in enumerate(atoms):
            am[i, i] = atom

        G = nx.from_numpy_matrix(am)

        return G

    def show_atom_number(self, mol, label):
        """
        Adds atom index labels as a property to each atom in an RDKit molecule,
        making them visible when drawing the molecule.

        This utility function is essential for visualizing the atom indices
        on the molecular diagram, which is crucial for interpreting the
        results of torsional angle detection.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The RDKit Mol object to which atom index labels will be added.

        label : str
            The **name of the atomic property** to set for displaying atom
            numbers (e.g., 'molAtomMapNumber' is commonly used by RDKit drawing functions).
            The value of this property will be the string representation of the atom's index.

        Return
        ------
        mol : rdkit.Chem.rdchem.Mol
            The **modified RDKit Mol object** with the atom map numbers set
            as an atomic property.
        """
        
        for atom in mol.GetAtoms():
            atom.SetProp(label, str(atom.GetIdx()))
        return mol

    def draw_molecule(self, mol, name_pic, atoms_idx):
        """
        Generates and saves a PNG image of the RDKit molecule, with specific
        atoms highlighted and atom numbers displayed.

        This function provides a visual representation of the molecule,
        emphasizing the atoms involved in the detected torsional angles.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The RDKit Mol object to be drawn.

        name_pic : str
            The **base name** for the output PNG file. The file will be saved
            as `<name_pic>.png` in the current working directory.

        atoms_idx : List[int] or List[Tuple[int, ...]]
            A **list of atom indices** to be highlighted in the generated
            image. This list can contain individual atom indices or tuples
            of atom indices (e.g., representing the four atoms of a dihedral).

        Return
        -------
        None
            This function does not return any value; it generates a PNG file
            as a side effect.
        """
        
        self.show_atom_number(mol, 'molAtomMapNumber')

        Draw.MolToFile(mol, str(name_pic) + ".png",
                       includeAtomNumbers=True, highlightAtoms=atoms_idx)

    def plot_trs_ang(self, name):
        """
        Plots and saves PNG images for each detected torsional angle in the
        molecule, highlighting the atoms involved in each angle.

        This method iterates through all identified torsional angles and generates
        a separate image for each, providing a detailed visualization of the
        conformational points of interest.

        Parameters
        ----------
        name : str
            The **base name** for the output PNG files. Each file will be named
            as `<name>_mol_<index>.png`, where `<index>` corresponds to the
            sequence of the torsional angle.

        Return
        -------
        None
            This function does not return any value; it generates multiple PNG
            files as a side effect. A message indicating the number of generated
            figures is printed to the console.
        """
        
        rdkit_mol = self.validate_mol(self.mol)
        lst_tor_angl = self.seek_angles()

        for idx, values in enumerate(lst_tor_angl):
            self.draw_molecule(rdkit_mol, str(name) + '_' + 'mol_' + str(idx), values)

        print("Torsional angle figures have been printed. Total: {}".format(len(lst_tor_angl)))
