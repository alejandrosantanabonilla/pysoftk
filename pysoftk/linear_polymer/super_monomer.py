"""
Super Monomer Module.

This module provides the `Sm` class, which handles the chemical linking
of two RDKit molecules via a designated placeholder atom. It includes
methods for full 3D monomer generation, as well as highly optimized 2D
topological merging for ultra-large polymer assemblies.
"""

from rdkit import Chem 
from rdkit.Chem import AllChem  
from rdkit import RDLogger

# Disable unnecessary RDKit warnings to keep the console output clean
RDLogger.DisableLog('rdApp.*') 

class Sm:
    """
    Class to create a new combined molecule (dimer/monomer) from two provided 
    RDKit Mol objects by linking them at specified placeholder atoms.
    

    This class supports standard SMARTS-based reactions with 3D embedding, 
    as well as ultra-fast direct 2D graph manipulation for large-scale 
    polymerization algorithms.

    Attributes
    ----------
    mol_1 : rdkit.Chem.rdchem.Mol
        The first RDKit molecule object, containing the placeholder atom.
    mol_2 : rdkit.Chem.rdchem.Mol
        The second RDKit molecule object, also containing the placeholder atom.
    atom : str
        The atomic symbol of the placeholder atom (e.g., 'Br', 'At', 'X') 
        indicating the connection points.
    """

    __slots__ = ['mol_1', 'mol_2', 'atom']

    def __init__(self, mol_1: Chem.Mol, mol_2: Chem.Mol, atom: str):
        """
        Initializes the Sm class with two molecules and a placeholder atom.

        Parameters
        ----------
        mol_1 : rdkit.Chem.rdchem.Mol
            The first reactant molecule.
        mol_2 : rdkit.Chem.rdchem.Mol
            The second reactant molecule.
        atom : str
            The atomic symbol used as the reactive placeholder.
        """
        self.mol_1 = mol_1
        self.mol_2 = mol_2
        self.atom = atom

    def constructor(self) -> Chem.Mol:
        """
        Combines the two input molecules using a SMARTS-defined chemical reaction.

        The reaction replaces the placeholder atoms and joins the molecules.
        It converts the resulting products to SMILES and back to a Mol object 
        to naturally sanitize and standardize the chemical graph.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The resulting combined molecule as an RDKit Mol object.

        Raises
        ------
        ValueError
            If the reaction fails to yield valid products or SMILES strings.
        """
        # Construct the SMARTS reaction string cleanly using an f-string
        react = f"[*:1][{self.atom}].[*:2][{self.atom}]>>[*:1][*:2]"
        
        # Create a reaction object from the SMARTS string
        rxn = AllChem.ReactionFromSmarts(react)
        results = rxn.RunReactants([self.mol_1, self.mol_2])

        m3_smiles_list = []
        if results: 
            for products_tuple in results:
                for product_mol in products_tuple:
                    m3_smiles_list.append(Chem.MolToSmiles(product_mol))
        else:
            raise ValueError(f"Reaction failed. Ensure '{self.atom}' exists in both monomers.")

        if not m3_smiles_list:
            raise ValueError("No product SMILES generated.")

        mol4 = Chem.MolFromSmiles(m3_smiles_list[0])
        if mol4:
            Chem.SanitizeMol(mol4)
        else:
            raise ValueError(f"Could not create molecule from SMILES: {m3_smiles_list[0]}")

        return mol4

    def _bond_order(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Assigns bond orders to a molecule based on its own template and 
        adds explicit hydrogen atoms.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            An RDKit molecule object requiring bond order assignment.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The molecule with assigned bond orders and explicit hydrogens.
        """
        new_mol = AllChem.AssignBondOrdersFromTemplate(mol, mol) 
        new_mol_h = Chem.AddHs(new_mol)
        return new_mol_h

    def build_topology_only(self) -> Chem.Mol:
        """
        Highly optimized method to connect molecules via direct 2D graph manipulation.
        

        Bypasses the RDKit reaction engine and SMILES conversion, allowing 
        instantaneous merging of large polymer chains. It strictly skips 3D 
        embedding and Hydrogen addition to maintain O(1) performance.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The cleanly connected 2D RDKit Mol object.
        """
        # 1. Combine into one container (disconnected)
        combined = Chem.CombineMols(self.mol_1, self.mol_2)
        rw_mol = Chem.RWMol(combined)
        
        # 2. Foolproof placeholder search using atom symbols
        matches = [a.GetIdx() for a in rw_mol.GetAtoms() if a.GetSymbol() == self.atom]
        
        # 3. Sort placeholders by their original fragment
        num_atoms_a = self.mol_1.GetNumAtoms()
        placeholders_a = [idx for idx in matches if idx < num_atoms_a]
        placeholders_b = [idx for idx in matches if idx >= num_atoms_a]
        
        # Fallback to the traditional constructor if placeholders are missing
        if not placeholders_a or not placeholders_b:
            return self.constructor()

        # 4. Connect A to B (Connects the last placeholder of A to the first of B)
        idx_p_a = placeholders_a[-1] 
        idx_p_b = placeholders_b[0]
        
        neighbor_a = rw_mol.GetAtomWithIdx(idx_p_a).GetNeighbors()[0]
        neighbor_b = rw_mol.GetAtomWithIdx(idx_p_b).GetNeighbors()[0]
        
        # 5. Add Bond and Remove Placeholders safely
        rw_mol.AddBond(neighbor_a.GetIdx(), neighbor_b.GetIdx(), Chem.BondType.SINGLE)
        rw_mol.RemoveAtom(max(idx_p_a, idx_p_b))
        rw_mol.RemoveAtom(min(idx_p_a, idx_p_b))
        
        # 6. Sanitize graph so RDKit updates its internal valency tables
        Chem.SanitizeMol(rw_mol)
        
        return rw_mol.GetMol()

    def mon_to_poly(self) -> Chem.Mol:
        """
        Prepares a single monomer unit by performing the chemical reaction, 
        adding Hydrogens, assigning bond orders, and embedding into 3D space.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The processed 3D RDKit Mol object, ready for legacy polymerization.
        """
        mol = self.constructor()
        new_mol_h = self._bond_order(mol)
        
        status = AllChem.EmbedMolecule(new_mol_h)
        if status == -1:
            AllChem.EmbedMolecule(new_mol_h, useRandomCoords=True)
            
        return new_mol_h

    def monomer(self) -> Chem.Mol:
        """
        Produces the final capped monomer unit.

        This method generates the 3D monomer and dynamically replaces any 
        unreacted terminal placeholder atoms with Hydrogen atoms.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The final, Hydrogen-capped 3D monomer object.
        """
        mol_intermediate = self.mon_to_poly()
        rw_mol = Chem.RWMol(mol_intermediate)
        
        # Identify remaining placeholder atoms
        atoms_to_modify = [
            a.GetIdx() for a in rw_mol.GetAtoms() if a.GetSymbol() == self.atom
        ]

        # Convert placeholders to Hydrogen (Atomic Number 1)
        for atom_idx in atoms_to_modify:
            rw_mol.GetAtomWithIdx(atom_idx).SetAtomicNum(1) 

        mol_after_replacement = rw_mol.GetMol()
        if mol_after_replacement:
            Chem.SanitizeMol(mol_after_replacement)

        # Re-apply bond orders and re-embed to account for geometric changes
        final_mol_processed = self._bond_order(mol_after_replacement)
        status = AllChem.EmbedMolecule(final_mol_processed)
        if status == -1:
            AllChem.EmbedMolecule(final_mol_processed, useRandomCoords=True)

        return final_mol_processed
