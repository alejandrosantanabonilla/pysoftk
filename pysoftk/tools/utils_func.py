import os
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def get_file_extension(path):
    """Extracts the filename and extension from a given path.

    Args:
        path (str): The full path to the file.

    Returns:
        tuple[str, str]: A tuple containing two strings:
            - The filename without the extension.
            - The file extension (including the dot, e.g., ".txt"). 
              If the path has no extension, this will be an empty string.
    """
    file_name, file_extension = os.path.splitext(path)
    return file_name, file_extension


def pattern_recon(pattern):
    """Identifies unique characters in a string pattern and sorts them.

    The sorting is case-insensitive alphabetical. For example, "aBcAb"
    would result in ['A', 'B', 'c'] (or ['a', 'b', 'c'] depending on
    which case is preserved by `set` for mixed-case duplicates, then sorted
    alphabetically based on their uppercase representation).

    Args:
        pattern (str): The input string pattern (e.g., "ABCA").

    Returns:
        list[str]: A list of unique characters from the pattern, sorted
                   alphabetically (case-insensitively, preserving the original case
                   of the first occurrence encountered by set for each character).
                   For "bBaAcC", output could be ['A', 'a', 'B', 'b', 'C', 'c'] if all are unique,
                   or ['A', 'B', 'C'] if 'a' is treated as 'A' by set (which it isn't).
                   Actually, `set` is case-sensitive. So "bBaAcC" -> {'b', 'B', 'a', 'A', 'c', 'C'}.
                   `sorted(..., key=str.upper)` then sorts these as A,a,B,b,C,c.
    """
    return sorted(list(set(pattern)), key=lambda c: c.upper())

def pattern_repl(pattern, tup_repl):
    """Replaces characters in a pattern string based on a list of replacement tuples.

    The function iterates through each tuple in `tup_repl`. Each tuple is
    expected to contain an old character (or substring) to be replaced and a new string
    (the replacement). After all replacements, the modified pattern string is
    split by the "+" delimiter, and any empty strings resulting from the split
    are removed.

    Args:
        pattern (str): The initial string pattern (e.g., "A+B+A").
        tup_repl (list[tuple[str, str]]): A list of tuples, where each tuple
            contains two strings: (substring_to_replace, replacement_string).
            Example: [('A', 'mol1'), ('B', 'mol2')]

    Returns:
        list[str]: A list of strings derived from the modified pattern after
                   replacements and splitting by "+". Empty strings are filtered out.
                   For pattern="A+B+A" and tup_repl=[('A','X'),('B','Y')],
                   it becomes "X+Y+X", then returns ['X', 'Y', 'X'].
    """
    for old_val, new_val in tup_repl:
        pattern = pattern.replace(old_val, new_val)

    valid_seq = pattern.split("+")
    return list(filter(None, valid_seq))


def pattern_mol_seq(mols, pattern):
    """Generates a sequence of molecule names based on a pattern string and a list of molecule names.

    It first identifies unique characters in the `pattern` (e.g., "A", "B" from "A+B+A").
    Then, it creates replacement tuples by pairing these unique characters (sorted)
    with the molecule names provided in the `mols` list, in their respective orders.
    Finally, it uses `pattern_repl` to substitute these characters in the
    `pattern` with their corresponding molecule names and splits the result.

    Example:
        mols = ["water", "ethanol"]
        pattern = "A+B+A"
        Unique chars in pattern (sorted alphabetically, case-insensitively): ['A', 'B']
        Replacement tuples: [('A', "water"), ('B', "ethanol")]
        Resulting sequence: ["water", "ethanol", "water"]

    Args:
        mols (list[str]): A list of molecule names. The order should correspond to the
                          alphabetically sorted unique characters from the `pattern`.
        pattern (str): A string defining the sequence pattern, using characters
                       as placeholders for molecules (e.g., "A+B+A").

    Returns:
        list[str]: A list of molecule names arranged according to the pattern.
                   Returns an empty list or partially resolved list if `mols`
                   and unique characters in `pattern` do not align as expected.
    """
    unq_val = pattern_recon(pattern) # e.g., ['A', 'B'] for "A+B+A"
    
    # The zip function will stop when the shortest iterable is exhausted.
    # If len(unq_val) != len(mols), the pairing might be incomplete.
    # Consider adding a check or warning if lengths do not match.
    # For example:
    # if len(unq_val) > len(mols):
    #     print(f"Warning: Pattern requires {len(unq_val)} unique molecule types, but only {len(mols)} provided.")
    # elif len(unq_val) < len(mols):
    #     print(f"Warning: More molecules provided ({len(mols)}) than unique types in pattern ({len(unq_val)}).")

    tup_repl = tuple(zip(unq_val, mols)) # e.g., (('A', 'mol1'), ('B', 'mol2'))

    return pattern_repl(pattern, tup_repl)


def count_plholder(mol, atom):
    """Counts the occurrences of a specific placeholder atom type in an RDKit molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object to inspect.
        atom (str): The atomic symbol of the placeholder atom to count (e.g., "Br", "*").

    Returns:
        int: The total number of placeholder atoms of the specified type found
             in the molecule.
    """
    num_br = 0
    # Iterate through each atom in the molecule
    for atm_obj in mol.GetAtoms(): # Renamed 'atoms' to 'atm_obj' for clarity
        # Check if the atom's symbol matches the specified placeholder symbol
        if atm_obj.GetSymbol() == str(atom):
            num_br += 1
    return num_br


def atom_neigh(mol, atom):
    """Identifies all placeholder atoms of a given type and their direct neighbors.

    For each atom in the molecule that matches the specified placeholder `atom` symbol,
    this function finds all its neighboring atoms. It returns a flat list of tuples,
    where each tuple contains the index of a placeholder atom and the index of one
    of its neighbors. If a placeholder has multiple neighbors, multiple tuples will
    be generated for that placeholder, one for each neighbor.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
        atom (str): The atomic symbol of the placeholder atom (e.g., "Br", "*").

    Returns:
        list[tuple[int, int]]: A list of (placeholder_atom_index, neighbor_atom_index)
                               tuples. The list is flattened.
                               Example: If P1 is a placeholder bonded to N1a and N1b, and
                               P2 is a placeholder bonded to N2a, the result could be
                               [(idx_P1, idx_N1a), (idx_P1, idx_N1b), (idx_P2, idx_N2a)].
                               The order depends on RDKit's internal atom iteration.
    """
    neigh = []
    # Using RWMol is not necessary if the molecule is not modified in this function.
    # Iterating directly on 'mol' is sufficient for read-only operations.
    for atm_obj in mol.GetAtoms(): # Renamed 'atoms' to 'atm_obj' for clarity
        if atm_obj.GetSymbol() == str(atom):
            # For each neighbor of the current placeholder atom
            for nbr in atm_obj.GetNeighbors():
                neigh.append((atm_obj.GetIdx(), nbr.GetIdx()))
    return neigh


def tuple_bonds(lst_atm_neigh):
    """Generates pairs of atom indices for creating new bonds from a list of placeholder-neighbor pairs.

    This function processes an ordered list of (placeholder_atom_index,
    neighbor_atom_index) tuples. It extracts the neighbor indices from the "internal"
    tuples in this list (i.e., excluding the first and the last tuple from `lst_atm_neigh`
    when collecting neighbor indices) and then pairs these collected neighbor indices
    sequentially. This is typically used to define new bonds that will connect
    monomers in a polymer chain.

    The function assumes `lst_atm_neigh` is appropriately ordered (e.g., by
    placeholder index along a chain) for the desired connectivity.

    Example:
        If `lst_atm_neigh` = `[(P0,N0), (P1,N1), (P2,N2), (P3,N3), (P4,N4)]`.
        - It considers internal elements for neighbor collection: `(P1,N1)`, `(P2,N2)`, `(P3,N3)`.
        - It extracts their neighbor indices: `bond_idx` becomes `[N1, N2, N3]`.
        - It then pairs them: `[(N1, N2)]`. The last element `N3` is unpaired if `bond_idx` has odd length.

    Args:
        lst_atm_neigh (list[tuple[int, int]]): An ordered list of tuples, where
            each tuple is (placeholder_atom_index, neighbor_atom_index).
            This list might be generated from `atom_neigh` and subsequently sorted/ordered.

    Returns:
        list[tuple[int, int]]: A list of tuples, where each tuple contains two
            neighbor_atom_indices that are intended to be bonded together.
            For example, `[(neighbor_idx1, neighbor_idx2), (neighbor_idx3, neighbor_idx4)]`.
            Returns an empty list if `bond_idx` has fewer than two elements.
    """
    bond_idx = []
    # The loop `range(1, len(lst_atm_neigh) - 1)` processes elements from index 1 up to index len-2.
    # This means it requires at least 3 elements in lst_atm_neigh for the loop to execute.
    # If len(lst_atm_neigh) < 3, bond_idx remains empty.
    for i in range(1, len(lst_atm_neigh) - 1):
        _placeholder_idx, neighbor_idx = lst_atm_neigh[i] # Original: a,b = lst_atm_neigh[i]
        bond_idx.append(neighbor_idx) # Original: bond_idx.append(b)

    # This pairs consecutive elements from bond_idx.
    # e.g., if bond_idx = [N1, N2, N3, N4], result is [(N1,N2), (N3,N4)].
    # If bond_idx = [N1, N2, N3], result is [(N1,N2)].
    return [(bond_idx[j], bond_idx[j+1])
            for j in range(0, len(bond_idx) - (len(bond_idx) % 2), 2)]


def create_pol(mol, atom, tpb):
    """Creates specified bonds in a molecule and removes internal placeholder atoms.

    This function modifies the input RDKit molecule `mol` by:
    1. Converting `mol` to an `RWMol` object to allow modifications.
    2. Adding new single bonds between atom pairs specified in `tpb`.
    3. Identifying all atoms matching the placeholder `atom` symbol (e.g., "Br", "*").
    4. Removing all occurrences of these placeholder atoms *except* for those
       that are effectively terminal (i.e., the ones with the smallest and
       largest atom indices among all identified placeholders, after sorting).

    This is typically used in polymer construction to form bonds between monomer
    units (connected via placeholders) and then remove the (now internal) placeholder
    atoms. Terminal placeholders might be retained.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object. This is often
            a single `Mol` object containing disconnected fragments or monomers
            linked by placeholders. It will be modified.
        atom (str): The atomic symbol of the placeholder atom (e.g., "Br", "*").
                    This symbol is used to identify atoms for removal.
        tpb (list[tuple[int, int]]): A list of tuples, where each tuple
            `(idx1, idx2)` specifies the atom indices to be connected by a
            new single bond.

    Returns:
        rdkit.Chem.rdchem.RWMol: The modified RDKit molecule object (as an RWMol)
                                 with new bonds added and specified internal
                                 placeholders removed.
    """
    rwmol1 = Chem.RWMol(mol) # Create an editable version of the molecule

    # Add specified bonds
    for ini, fin in tpb:
        rwmol1.AddBond(ini, fin, Chem.BondType.SINGLE)

    # Prepare query for placeholder atom
    # For generic wildcards like '*', MolFromSmarts is more robust.
    if atom == '*':
        plc_query = Chem.MolFromSmarts("[*]")
    else:
        # For specific elements like "Br", "Cl", MolFromSmiles is usually fine.
        plc_query = Chem.MolFromSmiles(atom)

    if plc_query is None:
        # If the placeholder symbol is invalid and a query molecule cannot be created.
        # Depending on desired behavior, could raise an error or return rwmol1 as is.
        # For now, printing a warning.
        print(f"Warning: Could not create query molecule from placeholder symbol: {atom}. No atoms will be removed.")
        return rwmol1 # Or raise ValueError(f"Invalid placeholder atom symbol: {atom}")

    # Find all occurrences of the placeholder atom
    # GetSubstructMatches returns a tuple of tuples; sum with an empty tuple flattens it.
    # Convert to list and sort by atom index.
    placeholder_indices = sorted(list(sum(rwmol1.GetSubstructMatches(plc_query), ())))

    # Remove "internal" placeholder atoms.
    # "Internal" are those not at the ends of the sorted list of placeholder indices.
    # This means if there are 2 or fewer placeholders, none are considered "internal" by this logic.
    if len(placeholder_indices) > 2:
        # Indices to remove are all placeholder indices except the first and the last one.
        # Atoms must be removed in reverse order of their indices to avoid re-indexing issues
        # during removal.
        indices_to_remove = sorted(placeholder_indices[1:-1], reverse=True)
        for i in indices_to_remove:
            rwmol1.RemoveAtom(i)
            
    return rwmol1
