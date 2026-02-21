import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Dict, Optional

# RDKit Imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # Suppress RDKit warnings

# Open Babel and Pysoftk Imports
from openbabel import pybel as pb
from pysoftk.tools.utils_ob import ff_ob_relaxation
from pysoftk.format_printers.format_mol import Fmt

class Monomer:
    """
    Base class for a monomer unit.

    Represents a single molecular building block, holding its RDKit
    molecule object and information about its connection points.
    """
    def __init__(self, id: int, generation: int, smiles: str, attachment_map_numbers: List[int]):
        """
        Initializes the Monomer.

        Args:
            id (int): A unique identifier for this monomer instance.
            generation (int): The generation layer this monomer belongs to (e.g., 0 for core).
            smiles (str): The SMILES string of the monomer, which must include
                          dummy atoms with map numbers (e.g., "[*:1]") for
                          each attachment point.
            attachment_map_numbers (List[int]): A list of the integer map numbers
                                                (e.g., [1, 2]) corresponding to the
                                                dummy atoms in the SMILES string.
        """
        self.id = id
        self.generation = generation
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"Invalid SMILES string for monomer {id}: {smiles}")

        # Add hydrogens and generate an initial 3D conformer for the monomer
        self.mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(self.mol, AllChem.ETKDGv3())

        self.attachment_map_numbers = attachment_map_numbers
        # Dictionary to store connections to other Monomer objects
        self.connections: Dict[int, 'Monomer'] = {}
        # The atom index in *this* monomer that connects to its parent
        self.parent_connection_atom_idx: Optional[int] = None
        # Maps local atom indices to global indices in the final combined molecule
        self.global_atom_indices: Dict[int, int] = {}
        # Stores the "real" atom indices that the dummy atoms are bonded to
        self.attachment_atom_indices = self._get_attachment_atom_indices()

    def _get_attachment_atom_indices(self) -> List[int]:
        """
        Finds the real atom indices that the dummy atoms are connected to.

        This helper function iterates through the RDKit molecule to find the
        dummy atoms (e.g., "[*:1]") and then identifies their neighboring
        atom, which is the actual point of connection.

        Returns:
            List[int]: A list of the RDKit atom indices for the connection points,
                       in the same order as `self.attachment_map_numbers`.
        """
        indices = {}
        for atom in self.mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetIntProp("molAtomMapNumber")
                if map_num in self.attachment_map_numbers:
                    # Assume dummy atom has only one neighbor
                    indices[map_num] = atom.GetNeighbors()[0].GetIdx()
        # Return indices in the order specified by the input list
        return [indices[mn] for mn in self.attachment_map_numbers if mn in indices]

class CustomMonomer(Monomer):
    """
    A specialized Monomer for dendrimer branching units.

    This class distinguishes between a single parent attachment point
    and multiple child attachment points, which is essential for
    building the directed graph of a dendrimer.
    """
    def __init__(self, id: int, generation: int, smiles: str,
                 parent_attachment_map_num: Optional[int],
                 child_attachment_map_nums: List[int]):
        """
        Initializes the CustomMonomer.

        Args:
            id (int): A unique identifier for this monomer instance.
            generation (int): The generation layer this monomer belongs to.
            smiles (str): The SMILES string of the branching unit.
            parent_attachment_map_num (Optional[int]): The map number for the
                                                       single connection to the
                                                       parent monomer. This is
                                                       `None` for the core.
            child_attachment_map_nums (List[int]): A list of map numbers for
                                                   the new attachment points
                                                   that will connect to the
                                                   next generation.
        """
        # Combine all map numbers for the base Monomer class
        all_maps = child_attachment_map_nums + ([parent_attachment_map_num] if parent_attachment_map_num else [])
        super().__init__(id, generation, smiles, all_maps)

        self.child_attachment_map_nums = child_attachment_map_nums
        # Maps child map numbers directly to their RDKit atom indices
        self.child_attachment_atom_indices_map: Dict[int, int] = {}

        # Create a quick lookup for all map numbers to their atom indices
        map_to_idx = dict(zip(self.attachment_map_numbers, self.attachment_atom_indices))

        if parent_attachment_map_num:
            self.parent_connection_atom_idx = map_to_idx.get(parent_attachment_map_num)

        # Populate the child-specific map
        for map_num in self.child_attachment_map_nums:
            self.child_attachment_atom_indices_map[map_num] = map_to_idx.get(map_num)

class SurfaceMonomer(Monomer):
    """
    A specialized Monomer for terminal surface (capping) groups.

    This monomer only has a single parent connection and no child
    connections, effectively terminating a branch.
    """
    def __init__(self, id: int, generation: int, smiles: str, parent_attachment_map_num: int):
        """
        Initializes the SurfaceMonomer.

        Args:
            id (int): A unique identifier for this monomer instance.
            generation (int): The generation layer this monomer belongs to.
            smiles (str): The SMILES string of the capping group.
            parent_attachment_map_num (int): The map number for the
                                             single connection to the
                                             parent monomer.
        """
        super().__init__(id, generation, smiles, [parent_attachment_map_num])
        # The first (and only) attachment index is the parent connection
        self.parent_connection_atom_idx = self.attachment_atom_indices[0]

class DendrimerBuilder:
    """
    A class to construct, optimize, and manage dendrimer molecules.

    This class orchestrates the entire workflow, from processing
    monomer definitions to building the connectivity graph, constructing
    the 2D molecule, and finally optimizing the 3D structure using
    external libraries (pysoftk/Open Babel).
    """

    def __init__(self, core_data: dict, branch_data: dict, generations: int, surface_data: Optional[dict] = None):
        """
        Initializes the DendrimerBuilder.

        Args:
            core_data (dict): A dictionary defining the core (G0) monomer.
                Expected keys: 'smiles', 'child_attachment_map_nums'.
            branch_data (dict): A dictionary defining the branching monomer.
                Expected keys: 'smiles', 'parent_attachment_map_num',
                'child_attachment_map_nums'.
            generations (int): The number of branching generations to build
                               (e.g., 2 for a G2 dendrimer).
            surface_data (Optional[dict]): An optional dictionary defining the
                surface capping group. If provided, the dendrimer will be
                capped. Expected keys: 'smiles', 'parent_attachment_map_num'.
        """
        self.core_data = core_data
        self.branch_data = branch_data
        self.generations = generations
        self.surface_data = surface_data
        
        # This list will hold the abstract graph of connected Monomer objects
        self.monomer_graph: List[Monomer] = []
        # This will hold the final, 3D-optimized molecule
        self.optimized_mol: Optional[pb.Molecule] = None

    def build(self, force_field: str = "UFF", steps: int = 1000, ff_thr: float = 0.1):
        """
        Orchestrates the full dendrimer construction and optimization process.

        This is the main method to run the workflow. It calls internal
        methods to generate the graph, build the 2D molecule, and
        run the 3D optimization.

        Args:
            force_field (str, optional): The Open Babel force field to use for
                                         optimization (e.g., "UFF", "MMFF94").
                                         Defaults to "UFF".
            steps (int, optional): The maximum number of optimization steps.
                                   Defaults to 1000.
            ff_thr (float, optional): The energy convergence threshold for the
                                      force field optimization. Defaults to 0.1.
        """
        print("Step 1: Generating abstract dendrimer graph...")
        self._generate_graph()
        
        print("Step 2: Building and optimizing 3D molecule...")
        # This internal method handles the RDKit build, conversion, and optimization
        self.optimized_mol = self._build_and_optimize_mol(force_field, steps, ff_thr)
        
        if self.optimized_mol:
            print("Dendrimer construction complete.")
        else:
            print("Error: Dendrimer construction failed.")

    def get_optimized_mol(self) -> Optional[pb.Molecule]:
        """
        Returns the final, optimized molecule as a Pybel object.

        Returns:
            Optional[pb.Molecule]: The 3D optimized Pybel molecule,
                                   or `None` if the build process failed.
        """
        return self.optimized_mol

    def plot_topology(self, output_filename: str = "dendrimer_topology.png"):
        """
        Plots the abstract topological graph of the dendrimer using NetworkX.

        This visualization shows the monomer connectivity and generation,
        not the chemical structure.

        Args:
            output_filename (str, optional): The file path to save the
                                             plot image.
                                             Defaults to "dendrimer_topology.png".
        """
        if not self.monomer_graph:
            print("Cannot plot topology. Generate the graph first by running .build()")
            return
            
        graph = nx.Graph()
        labels = {m.id: f"G{m.generation}" for m in self.monomer_graph}
        for m in self.monomer_graph:
            graph.add_node(m.id)
            for child in m.connections.values():
                graph.add_edge(m.id, child.id)

        pos = nx.kamada_kawai_layout(graph)
        node_colors = [m.generation for m in self.monomer_graph]

        fig, ax = plt.subplots(figsize=(10, 8))
        nx.draw(graph, pos, ax=ax, labels=labels, node_size=700,
                node_color=node_colors, cmap=plt.cm.viridis, font_size=8)
        
        ax.set_title("Dendrimer Topological Structure")
        plt.savefig(output_filename)
        print(f"Topology plot saved to '{output_filename}'")
        plt.close(fig)

    def _generate_graph(self):
        """
        (Internal) Creates the abstract graph of connected Monomer objects.

        This method builds the dendrimer layer by layer in memory,
        populating `self.monomer_graph` with the connected
        Monomer, CustomMonomer, and SurfaceMonomer objects.
        """
        monomer_count = 0
        # Create the G0 core
        core = CustomMonomer(0, 0, self.core_data['smiles'], None, self.core_data['child_attachment_map_nums'])
        self.monomer_graph = [core]
        
        # Build branching generations
        for gen in range(1, self.generations + 1):
            parents = [m for m in self.monomer_graph if m.generation == gen - 1]
            for parent in parents:
                if isinstance(parent, CustomMonomer):
                    # Add a new child for each of the parent's child attachment points
                    for map_num in parent.child_attachment_map_nums:
                        monomer_count += 1
                        child = CustomMonomer(monomer_count, gen, **self.branch_data)
                        parent.connections[map_num] = child
                        self.monomer_graph.append(child)
        
        # Add surface capping groups if specified
        if self.surface_data:
            gen_to_cap = self.generations
            surface_parents = [m for m in self.monomer_graph if m.generation == gen_to_cap]
            for parent in surface_parents:
                if isinstance(parent, CustomMonomer):
                    # Add a new cap for each of the final generation's attachment points
                    for map_num in parent.child_attachment_map_nums:
                        monomer_count += 1
                        cap = SurfaceMonomer(monomer_count, parent.generation + 1, **self.surface_data)
                        parent.connections[map_num] = cap
                        self.monomer_graph.append(cap)
    
    def _build_and_optimize_mol(self, force_field: str, steps: int, ff_thr: float) -> Optional[pb.Molecule]:
        """
        (Internal) Handles the full molecule construction and optimization.

        This method:
        1. Builds the complete 2D molecular graph in RDKit (RWMol).
        2. Forms all inter-monomer bonds.
        3. Finalizes and sanitizes the RDKit molecule.
        4. Converts the RDKit molecule to a Pybel molecule.
        5. Generates an initial 3D structure.
        6. Runs the pysoftk force field optimization.

        Args:
            force_field (str): The Open Babel force field to use.
            steps (int): The maximum number of optimization steps.
            ff_thr (float): The energy convergence threshold.

        Returns:
            Optional[pb.Molecule]: The optimized Pybel molecule, or `None` on failure.
        """
        rw_mol = Chem.RWMol()
        
        # 1. Add all atoms and internal bonds from each monomer to the RWMol
        for monomer in self.monomer_graph:
            mol_copy = Chem.Mol(monomer.mol)
            for atom in mol_copy.GetAtoms():
                if not atom.HasProp("molAtomMapNumber"):
                    # Add non-dummy atoms and store their new global index
                    new_idx = rw_mol.AddAtom(atom)
                    monomer.global_atom_indices[atom.GetIdx()] = new_idx
            for bond in mol_copy.GetBonds():
                b, e = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                # Add internal bonds only between atoms that were added
                if b in monomer.global_atom_indices and e in monomer.global_atom_indices:
                    rw_mol.AddBond(monomer.global_atom_indices[b], monomer.global_atom_indices[e], bond.GetBondType())

        # 2. Add all inter-monomer bonds
        for parent in self.monomer_graph:
            if isinstance(parent, CustomMonomer):
                for map_num, child in parent.connections.items():
                    # Find the global indices for the parent's and child's connection atoms
                    p_idx = parent.global_atom_indices[parent.child_attachment_atom_indices_map[map_num]]
                    c_idx = child.global_atom_indices[child.parent_connection_atom_idx]
                    rw_mol.AddBond(p_idx, c_idx, BondType.SINGLE)
        
        # 3. Finalize RDKit molecule
        try:
            final_mol_rdkit = rw_mol.GetMol()
            Chem.SanitizeMol(final_mol_rdkit)
        except Exception as e:
            print(f"  Error during RDKit molecule sanitization: {e}")
            return None
        
        # 4. Convert to Pybel object
        try:
            mol_block = Chem.MolToMolBlock(final_mol_rdkit)
            mol_pybel = pb.readstring('mol', mol_block)
        except Exception as e:
            print(f"  Error converting RDKit molecule to Pybel: {e}")
            return None

        # 5. Generate 3D structure and add hydrogens
        mol_pybel.addh()
        mol_pybel.make3D()
        
        # 6. Optimize using pysoftk's ff_ob_relaxation function
        print(f"    Optimizing with pysoftk using {force_field} force field ({steps} steps, threshold {ff_thr})...")
        try:
            # The function ff_ob_relaxation takes steps as the 3rd argument (relax_iterations)
            optimized_mol_pybel = ff_ob_relaxation(mol_pybel, force_field, steps, ff_thr)
            print("    Optimization complete.")
        except Exception as e:
            print(f"    An exception occurred during pysoftk optimization: {e}")
            return None
            
        return optimized_mol_pybel
