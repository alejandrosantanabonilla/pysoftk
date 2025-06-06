from rdkit import Chem 
from rdkit.Chem import AllChem  
from rdkit.Chem import rdDistGeom as molDG 

# Disable unnecessary RDKit warnings to keep the output clean
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # Supresses common RDKit informational messages

class Sm(object):
    """
    Class to create a new combined molecule (monomer) from two provided RDKit Mol objects
    by linking them at specified placeholder atoms.

    Example
    -------
    # Assuming mol1 and mol2 are RDKit Mol objects and 'X' is the placeholder atom
    # from rdkit import Chem
    # mol1 = Chem.MolFromSmiles('Cc1ccccc1[X]') # Example molecule 1 with placeholder X
    # mol2 = Chem.MolFromSmiles('O=C(O)c1ccc(N)cc1[X]') # Example molecule 2 with placeholder X
    # combiner = Sm(mol1, mol2, 'X')
    # new_monomer = combiner.monomer()
    # print(Chem.MolToSmiles(new_monomer))

    Note
    ----
    This class requires RDKit to be installed.
    The placeholder atom should be a unique atom type (e.g., a dummy atom like '*')
    or an atom that is intended to be replaced during the reaction.
    """

    # __slots__ can provide memory optimization by pre-declaring instance attributes.
    # This can be beneficial if many instances of this class are created.
    __slots__ = ['mol_1', 'mol_2' , 'atom']

    def __init__(self, mol_1, mol_2, atom):
        """
        Initializes the Sm class with two molecules and a placeholder atom.

        Parameters
        ----------
        mol_1 : rdkit.Chem.rdchem.Mol
            The first RDKit molecule object. This molecule should contain the placeholder atom.
        mol_2 : rdkit.Chem.rdchem.Mol
            The second RDKit molecule object. This molecule should also contain the placeholder atom.
        atom : str
            The atomic symbol of the placeholder atom (e.g., 'X', 'R', '*')
            which indicates the connection points in mol_1 and mol_2.
        """
        
        # Store the input molecules and placeholder atom as instance attributes
        self.mol_1 = mol_1
        self.mol_2 = mol_2
        self.atom = atom

    def constructor(self):
        """
        Combines the two input molecules (mol_1 and mol_2) using a chemical reaction
        defined by a SMARTS pattern. The reaction replaces the placeholder atoms
        and joins the molecules.

        The SMARTS reaction pattern '[*:1][AtomSymbol].[*:2][AtomSymbol]>>[*:1][*:2]'
        essentially means:
        - Find a fragment with any atom labeled as map number 1 ([*:1]) connected to the placeholder atom ([AtomSymbol]).
        - Find another fragment with any atom labeled as map number 2 ([*:2]) connected to the placeholder atom ([AtomSymbol]).
        - Combine these two fragments by forming a bond between the atom mapped as 1 and the atom mapped as 2,
          effectively removing the placeholder atoms.

        Returns
        -------
        mol4 : rdkit.Chem.rdchem.Mol
            The resulting combined molecule as an RDKit Mol object.
            Returns the first product of the reaction.
        """
        
        # Retrieve instance attributes for local use
        mol_1, mol_2, atom = self.mol_1, self.mol_2, self.atom

        # Construct the SMARTS reaction string.
        # [*:1] and [*:2] are atom maps used to define which atoms form the new bond.
        # The placeholder atom (self.atom) is specified as the atom to be removed/replaced.
        # For example, if atom is 'X', the pattern becomes '[*:1][X].[*:2][X]>>[*:1][*:2]'
        # This means: find an X atom attached to something (map 1), find another X atom
        # attached to something else (map 2), and then connect map 1 to map 2, removing the X's.
        patt = ['[*:1]', '[', str(atom), ']', '.', '[*:2]',
                '[', str(atom), ']', '>>', '[*:1]', '[*:2]']
        react = "".join(patt)  # Join the list into a single string

        # Create a reaction object from the SMARTS string
        rxn = AllChem.ReactionFromSmarts(str(react))

        # Run the reaction with mol_1 and mol_2 as reactants
        # RunReactants returns a list of lists of product molecules.
        # Each inner list represents a set of products from one reaction event.
        results = rxn.RunReactants([mol_1, mol_2])

        # Initialize a list to store SMILES strings of the products
        m3_smiles_list = []
        # Iterate through the reaction results
        # `results` is a tuple of tuples of product molecules.
        # Example: (((product1_mol, product2_mol), (product1_mol_variant2, ...)), ...)
        # We are typically interested in the first set of products from the first reaction outcome.
        if results: # Check if any products were formed
            for products_tuple in results: # Each products_tuple is a tuple of Mol objects
                for product_mol in products_tuple:
                    # Convert each product molecule to its SMILES representation
                    m3_smiles_list.append(Chem.MolToSmiles(product_mol))
        else:
            # Handle the case where the reaction yields no products
            # This could be due to various reasons, e.g., reactants don't match the pattern.
            raise ValueError("Reaction did not yield any products. Check input molecules and placeholder atom.")

        # Check if any product SMILES were generated
        if not m3_smiles_list:
            raise ValueError("No product SMILES generated. The reaction might have failed or produced no valid molecules.")

        # Create a new RDKit Mol object from the SMILES string of the first product
        # This assumes that the desired product is the first one generated.
        # For more complex reactions, more sophisticated product selection might be needed.
        mol4 = Chem.MolFromSmiles(m3_smiles_list[0])

        # It's good practice to sanitize the molecule to ensure correct chemistry
        if mol4:
            Chem.SanitizeMol(mol4)
        else:
            raise ValueError(f"Could not create molecule from SMILES: {m3_smiles_list[0]}")


        return mol4

    def _bond_order(self, mol):
        """
        Assigns bond orders to a molecule based on a template (itself in this case)
        and then adds explicit hydrogen atoms. This is often used to ensure that
        the molecule has a chemically sensible representation, especially after
        manipulations like reactions.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            An RDKit molecule object whose bond orders need to be assigned and hydrogens added.

        Returns
        -------
        newMol_H : rdkit.Chem.rdchem.Mol
            The molecule with assigned bond orders and added explicit hydrogens.
        """
        
        # Assign bond orders using the molecule itself as a template.
        # This can help kekulize aromatic systems and infer bond orders if they are ambiguous.
        # Note: If `mol` comes from a reaction, its bond orders might already be set.
        # This step might be more relevant if `mol` was constructed from, e.g., an XYZ file.
        newMol = AllChem.AssignBondOrdersFromTemplate(mol, mol) # Using mol as its own template

        # Add explicit hydrogen atoms to satisfy valencies.
        # The `addCoords=True` option (default is False) could be considered if 3D coordinates are important
        # and need to be generated for the new hydrogens, though EmbedMolecule later handles this.
        newMol_H = Chem.AddHs(newMol)

        return newMol_H

    def mon_to_poly(self):
        """
        Prepares a single monomer unit for potential polymerization or further processing.
        It constructs the combined molecule, assigns bond orders, adds hydrogens,
        and generates a 3D conformation.

        This method essentially chains `constructor`, `_bond_order`, and `EmbedMolecule`.

        Returns
        -------
        newMol_H : rdkit.Chem.rdchem.Mol
            The processed RDKit Mol object, ready for use as a monomer,
            with 3D coordinates.
        """
        
        # Step 1: Combine the initial two molecules to form the basic monomer structure.
        mol = self.constructor()

        # Step 2: Assign bond orders and add explicit hydrogen atoms.
        newMol_H = self._bond_order(mol)

        # Step 3: Generate a 3D conformation for the molecule.
        # AllChem.EmbedMolecule attempts to generate a reasonable 3D structure.
        # It uses distance geometry methods.
        # `useRandomCoords=True` can sometimes help find better conformations if the default fails.
        # `maxAttempts` can be increased if embedding is difficult.
        # The function returns an status code (0 for success). It's good practice to check this.
        status = AllChem.EmbedMolecule(newMol_H)
        if status == -1:
            print(f"Warning: Could not embed molecule: {Chem.MolToSmiles(newMol_H)}. Conformation may be poor or missing.")
            # Optionally, try with random coordinates if the first attempt fails
            status = AllChem.EmbedMolecule(newMol_H, useRandomCoords=True)
            if status == -1:
                print(f"Warning: Embedding with random coordinates also failed for {Chem.MolToSmiles(newMol_H)}.")

        return newMol_H

    def monomer(self):
        """
        Produces the final monomer unit.
        This function takes the molecule prepared by `mon_to_poly`,
        replaces any remaining placeholder atoms (specified by `self.atom`) with hydrogen atoms,
        and then re-processes (bond order assignment, H addition, embedding) the molecule.
        This is useful if the placeholder atom was not fully consumed or if it needs to be
        explicitly capped with hydrogen in the final monomer.

        Returns
        -------
        final_mol : rdkit.Chem.rdchem.Mol
            The resulting final monomer as an RDKit Mol object with 3D coordinates.
        """
        
        # Get the placeholder atom symbol
        placeholder_symbol = str(self.atom)

        # Step 1: Generate the initial combined and 3D embedded molecule.
        # This molecule might still contain instances of the placeholder atom if the
        # reaction in `constructor` was designed to leave some, or if it's a different
        # placeholder atom intended for this step.
        mol_intermediate = self.mon_to_poly()

        # Step 2: Create an editable version of the molecule (RWMol).
        # RWMol allows addition and removal of atoms and bonds.
        rw_mol = Chem.RWMol(mol_intermediate)

        # Step 3: Iterate through all atoms in the molecule.
        # If an atom matches the placeholder symbol, change its atomic number to 1 (Hydrogen).
        atoms_to_modify = [] # Store indices to avoid issues with modifying during iteration
        for atom_obj in rw_mol.GetAtoms():
            if atom_obj.GetSymbol() == placeholder_symbol:
                atoms_to_modify.append(atom_obj.GetIdx())

        for atom_idx in atoms_to_modify:
            rw_mol.GetAtomWithIdx(atom_idx).SetAtomicNum(1) # Change to Hydrogen

        # Step 4: Convert the RWMol back to a regular Mol object.
        # SanitizeMol is crucial here to update valencies, aromaticity, etc., after changing atom types.
        mol_after_replacement = rw_mol.GetMol()
        if mol_after_replacement:
            try:
                Chem.SanitizeMol(mol_after_replacement)
            except Exception as e:
                print(f"Error during sanitization after replacing placeholder with H: {e}")
                # Potentially return mol_intermediate or raise an error, depending on desired behavior
                return mol_intermediate # Or handle error appropriately

        # Step 5: Re-apply bond order assignment and add hydrogens.
        # This is important because changing atom types (e.g., placeholder to H)
        # can affect valency and bonding.
        final_mol_processed = self._bond_order(mol_after_replacement)

        # Step 6: Re-embed the molecule to get updated 3D coordinates.
        # The geometry might change significantly after replacing atoms and adding hydrogens.
        status = AllChem.EmbedMolecule(final_mol_processed)
        if status == -1:
            print(f"Warning: Could not re-embed molecule after H replacement: {Chem.MolToSmiles(final_mol_processed)}. Conformation may be poor or missing.")
            status = AllChem.EmbedMolecule(final_mol_processed, useRandomCoords=True)
            if status == -1:
                print(f"Warning: Re-embedding with random coords also failed for {Chem.MolToSmiles(final_mol_processed)}.")

        return final_mol_processed
