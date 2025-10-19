import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import copy
from collections import OrderedDict

from pysoftk.linear_polymer import linear_polymer as lp
from pysoftk.linear_polymer import super_monomer as sm
from pysoftk.topologies.ring import *

from pysoftk.tools.utils_rdkit import *
from pysoftk.tools.utils_func import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

class Db:
    """
    A class specifically designed for the **computational construction of diblock copolymers**
    from two distinct RDKit molecular monomers (Block A and Block B).

    This class provides a streamlined workflow to:
    1. Generate individual polymer blocks (Block A and Block B) of specified lengths.
    2. Concatenate these two blocks to form a single diblock copolymer structure.
    3. Perform subsequent 3D embedding and force field optimization to obtain a stable
       and chemically sensible molecular geometry.

    Examples
    ---------
    >>> from rdkit import Chem
    >>> from pysoftk.linear_polymer.linear_polymer import Lp
    >>> from pysoftk.block_copolymers.diblock import Db
    >>> # Define monomer A and monomer B with placeholder atoms (e.g., [Pt])
    >>> monomer_A = Chem.MolFromSmiles("C([Pt])C")
    >>> monomer_B = Chem.MolFromSmiles("CC([Pt])C")
    >>>
    >>> # Initialize the Db class with the two monomers and the placeholder atom
    >>> diblock_builder = Db(ma=monomer_A, mb=monomer_B, atom="Pt")
    >>>
    >>> # Create a diblock copolymer with 5 units of Block A and 3 units of Block B
    >>> final_diblock_mol = diblock_builder.diblock_copolymer(
    ...     len_block_A=5, len_block_B=3, force_field="MMFF"
    ... )
    >>> print(f"Diblock copolymer atoms: {final_diblock_mol.GetNumAtoms()}")

    Note
    -----
    The **RDKit package must be installed** in your Python environment for this
    class to function correctly. Ensure that the placeholder atom specified (e.g., "Pt")
    is consistently used in both monomer A and monomer B. This class relies on
    functionalities from `pysoftk.linear_polymer` and `pysoftk.tools`.
    """

    __slots__ = ['ma', 'mb', 'atom']

    def __init__(self, ma, mb, atom):
        """
        Initializes the `Db` class by defining the two distinct monomer RDKit
        molecules (monomer A and monomer B) and the common placeholder atom
        used for linking them.

        These input molecules serve as the fundamental repeating units for
        the respective blocks within the diblock copolymer.

        Parameters
        -----------
        ma : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing **monomer A**, which will
            form the first block of the diblock copolymer. This molecule
            must contain a designated placeholder atom.

        mb : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing **monomer B**, which will
            form the second block of the diblock copolymer. This molecule
            must also contain a designated placeholder atom.

        atom : str
            The **atomic symbol** (e.g., "Pt", "Au") of the placeholder atom.
            This specific atom type must be present in both `ma` and `mb`
            and will be used to establish bonds during the polymerization
            process, connecting individual monomers within each block and
            then connecting Block A to Block B.
        """

        self.ma = ma
        self.mb = mb
        self.atom = atom

    def diblock_copolymer(self, len_block_A, len_block_B,
                          force_field="MMFF", relax_iterations=100, rot_steps=1):
        """
        Constructs a complete **diblock copolymer** by first synthesizing
        Block A and Block B independently to specified lengths, and then
        covalently linking them. The final copolymer undergoes 3D embedding
        and force field optimization to generate a stable conformation.

        This function orchestrates the entire diblock copolymer synthesis process,
        from individual block formation to the final optimized structure.

        Parameters
        -----------
        len_block_A : int
            The **desired number of monomer A units** that will comprise
            the first block (Block A) of the diblock copolymer.

        len_block_B : int
            The **desired number of monomer B units** that will comprise
            the second block (Block B) of the diblock copolymer.

        force_field : str, optional
            The **name of the force field** to apply during the geometry
            optimization steps for both individual blocks and the final
            diblock copolymer. Common choices include "MMFF" (which defaults
            to MMFF94) or "UFF". Defaults to "MMFF".

        relax_iterations : int, optional
            The **maximum number of iterations** for the force field relaxation
            algorithm. A higher value generally leads to better energy minimization
            but increases computation time. This value is applied to the
            optimization of both blocks and the final copolymer. Defaults to 100.

        rot_steps : int, optional
            The **number of rotational steps** to perform during conformational
            searching. This helps in exploring different low-energy conformers
            of the polymer. Defaults to 1.

        Return
        -------
        diblock : openbabel.pybel.Molecule
            An **OpenBabel Pybel Molecule object** representing the
            fully constructed and geometry-optimized diblock copolymer.
            This object contains the 3D coordinates and bonding information
            of the final polymer.
        """

        ma = self.ma
        mb = self.mb
        atom = self.atom

        string_1 = 'A' * len_block_A
        monomer = Pt(string_1, [ma], str(atom)).pattern_block_poly(relax_iterations, force_field=str(force_field), swap_H=False)


        string_2 = 'B' * len_block_B
        monomer2 = Pt(string_2, [mb], str(atom)).pattern_block_poly(relax_iterations, force_field=str(force_field), swap_H=False)

        mon_temp = monomer.write("mol")
        mon2_temp = monomer2.write("mol")

        mon = Chem.MolFromMolBlock(mon_temp)
        mon2 = Chem.MolFromMolBlock(mon2_temp)

        diblock = sm.Sm(mon, mon2, str(atom)).monomer()

        last_mol = check_proto(diblock, force_field, relax_iterations, rot_steps)

        return last_mol

class Pt:
    """
    A versatile class for creating **patterned polymers** based on a
    user-defined sequence of RDKit molecular building blocks.

    This class enables the construction of complex copolymers by
    allowing a specific sequence (pattern) of different monomers to be
    arranged into a linear chain. It supports 3D embedding and force field
    optimization of the resulting polymer.

    Examples
    ---------
    >>> from rdkit import Chem
    >>> from pysoftk.block_copolymers.diblock import Pt
    >>> # Define two distinct monomers (e.g., A and B) with a placeholder atom
    >>> monomer_A = Chem.MolFromSmiles("C([Pt])C")
    >>> monomer_B = Chem.MolFromSmiles("N([Pt])N")
    >>>
    >>> # Define a repeating pattern
    >>> pattern = "ABABAB"
    >>> # List of molecules corresponding to the pattern's unique characters
    >>> molecules = [monomer_A, monomer_B]
    >>>
    >>> # Initialize the Pt class
    >>> patterned_poly_builder = Pt(pattern=pattern, mols=molecules, atom="Pt")
    >>>
    >>> # Build the patterned polymer
    >>> final_patterned_mol = patterned_poly_builder.pattern_block_poly(
    ...     relax_iterations=200, force_field="UFF", swap_H=True
    ... )
    >>> print(f"Patterned polymer atoms: {final_patterned_mol.GetNumAtoms()}")

    Note
    -----
    The **RDKit package must be installed** in your Python environment.
    The `pattern` string should consist of characters that correspond to
    the order of molecules in the `mols` list (e.g., 'A' maps to `mols[0]`,
    'B' maps to `mols[1]`, etc.). Ensure consistency in the `atom` placeholder.
    """

    __slots__ = ['pattern', 'mols', 'atom']

    def __init__(self, pattern, mols, atom):
        """
        Initializes the `Pt` class with the desired polymerization pattern,
        the list of RDKit molecular building blocks, and the placeholder atom.

        Parameters
        -----------
        pattern : str
            A **string representing the sequence** in which the molecular
            building blocks should be arranged to form the polymer. Each
            character in the string corresponds to a specific molecule in
            the `mols` list (e.g., "ABAB" or "AAABBB").

        mols : list[rdkit.Chem.rdchem.Mol]
            A **list of RDKit Mol objects**, where each object represents a
            unique monomer type. The order of molecules in this list should
            correspond to the alphabetical order of characters used in the `pattern`
            string (e.g., `mols[0]` for 'A', `mols[1]` for 'B', and so on).
            Each monomer must contain the specified `atom` as a placeholder.

        atom : str
            The **atomic symbol** (e.g., "Au", "Sn") of the placeholder atom
            present in all the molecular objects within the `mols` list. This
            atom facilitates the merging and bonding of the different monomer
            units according to the `pattern`.
        """

        self.pattern = pattern
        self.mols = mols
        self.atom = atom

    def pattern_block_poly(self, relax_iterations=100, force_field="MMFF", swap_H=True, rot_steps=1, more_iter=10):
        """
        Constructs a **linear polymer according to a defined alphabetical pattern**
        of molecular building blocks, followed by geometry optimization.

        This function assembles the polymer by sequentially combining the
        specified monomers based on the input `pattern` string. It handles
        the creation of bonds between units, and then performs comprehensive
        3D structure optimization using a force field, with options for
        replacing placeholder atoms and further refinement.

        Parameters
        -----------
        relax_iterations : int, optional
            The **maximum number of iterations** for the force field relaxation
            procedure. A higher number typically leads to better convergence
            of the molecular geometry. Defaults to 100.

        force_field : str, optional
            The **name of the force field** to be used for molecular mechanics
            optimization. Supported options generally include "MMFF" (which maps
            to MMFF94) or "UFF". Defaults to "MMFF".

        swap_H : bool, optional
            A flag indicating whether the placeholder atoms (specified by `self.atom`)
            at the termini of the constructed polymer should be **replaced with
            hydrogen atoms**.
            * If `True` (default), placeholder atoms are converted to hydrogens
              before the final optimization.
            * If `False`, the placeholder atoms remain in the structure.

        rot_steps : int, optional
            The **number of rotational steps** to explore during conformational
            optimization. This parameter helps in finding different low-energy
            arrangements of the polymer chain by rotating around rotatable bonds.
            Defaults to 1.

        more_iter : int, optional
            A **multiplier for additional optimization iterations**. The
            `relax_iterations` value will be multiplied by this number to perform
            a more extensive final geometry optimization, potentially leading
            to a more stable conformation. Defaults to 10.

        Return
        --------
        mol : openbabel.pybel.Molecule
            A **PySoftK molecular object** (which wraps an OpenBabel Pybel Molecule)
            representing the fully constructed and geometry-optimized patterned polymer.
            This object contains the 3D coordinates and bonding information of the final structure.
        """
        
        pattern = self.pattern
        mols = self.mols
        atom = self.atom

        names = ['mol_{}+'.format(i) for i in range(1, len(mols)+1)]
        seq = pattern_mol_seq(names, pattern)

        numbers = ['mol_{}'.format(i) for i in range(1, len(mols)+1)]
        od_mols = OrderedDict(list(zip(numbers, mols)))

        list_mol = [od_mols[i] for i in seq]
        outmol = Chem.CombineMols(list_mol[0], list_mol[1])

        for i in range(2, len(list_mol)):
            outmol = Chem.CombineMols(outmol, list_mol[i])

        lst_ngh = atom_neigh(outmol, str(atom))
        tpb = tuple_bonds(lst_ngh)
        proto_pol = create_pol(outmol, str(atom), tpb)

        if swap_H:
            newMol_H = swap_hyd(proto_pol, relax_iterations, str(atom), force_field)

        if not swap_H:
            newMol_H = no_swap(proto_pol, relax_iterations, force_field)

        last_mol = check_proto(newMol_H, force_field, relax_iterations*more_iter, rot_steps)

        return last_mol
