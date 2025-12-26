from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom, rdMolTransforms
from rdkit.Chem import rdDistGeom as molDG

import numpy as np
from pysoftk.tools.utils_func import *
from pysoftk.tools.utils_rdkit import *

class Bd:
    """
    A class designed for the **automated construction of branched polymers**
    from user-defined molecular building blocks (a core and arm molecules).

    This class provides functionalities to attach multiple 'arm' molecules to a
    'core' molecule, forming a branched structure. It handles the bonding
    using a placeholder atom strategy and includes options for 3D embedding
    and force field optimization of the final branched polymer.

    Examples
    ---------
    >>> from rdkit import Chem
    >>> from pysoftk.branched_polymers.branched_polymer import Bd
    >>> # Define a core molecule with placeholder atoms (e.g., [Au])
    >>> core_smiles = "C([Au])C([Au])C"
    >>> core_mol = Chem.MolFromSmiles(core_smiles)
    >>> # Define an arm molecule with a placeholder atom
    >>> arm_smiles = "CC([Au])"
    >>> arm_mol = Chem.MolFromSmiles(arm_smiles)
    >>> # Initialize the Bd class
    >>> branched_poly_builder = Bd(core=core_mol, arm=arm_mol, atom="Au")
    >>> # Build the branched polymer
    >>> polymer_mol = branched_poly_builder.branched_polymer(force_field="MMFF", relax_iterations=200)
    >>> print(f"Branched polymer atoms: {polymer_mol.GetNumAtoms()}")
    >>> print(f"Branched polymer bonds: {polymer_mol.GetNumBonds()}")

    Note
    -----
    The **RDKit package must be installed** in your Python environment for this
    class to function. Ensure that the placeholder atom specified (e.g., "Au")
    is consistent across your core and arm molecules and is a valid RDKit atom.
    """

    __slots__ = ['core', 'arm', 'atom']

    def __init__(self, core, arm, atom):
        """
        Initializes the `Bd` class by defining the core molecule, the arm
        molecule, and the placeholder atom used for connecting them.

        These input molecules form the fundamental building blocks for the
        branched polymer structure.

        Parameters
        ----------
        core : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the central core molecule
            of the branched polymer. This molecule should contain one or more
            placeholder atoms to which the arm molecules will be attached.

        arm : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the arm molecule that will
            be repeatedly attached to the core. This molecule must also contain
            a placeholder atom for connection.

        atom : str
            The **atomic symbol** (e.g., "Au", "Pt") of the placeholder atom
            present in both the `core` and `arm` molecules. This atom will be
            used to identify attachment points and facilitate the formation of
            new bonds between the core and the arms.
        """

        self.core = core
        self.arm = arm
        self.atom = atom

    def merge_arms(self, core, arm, atom):
        """
        Attaches a single `arm` molecule to a `core` molecule by identifying
        and forming a bond between their respective placeholder atoms, then
        removing these placeholders.

        This function performs the fundamental operation of merging a branch
        onto the growing polymer structure. It combines the molecules, embeds
        them in 3D space, identifies the placeholder atoms, creates a new bond
        between their connected atoms, and finally removes the placeholder atoms.

        Parameters
        ----------
        core : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the current core structure
            (which could be the initial core or a partially built branched polymer)
            to which an arm molecule will be attached. It must contain at least
            one placeholder atom.

        arm : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the arm molecule to be attached.
            It must contain a placeholder atom.

        atom : str
            The **atomic symbol** of the placeholder atom used in both the `core`
            and `arm` molecules that will be involved in the bonding.

        Return
        -------
        mol : rdkit.Chem.rdchem.Mol
            A **new RDKit Mol object** representing the combined structure
            after the `arm` molecule has been successfully attached to the `core`,
            with the placeholder atoms removed and a new bond formed.
        """

        outmol = Chem.CombineMols(core, arm)
        outmol.UpdatePropertyCache(strict=False)
        AllChem.EmbedMolecule(outmol)

        new_bond = atom_neigh(outmol, str(atom))

        br_1, c_1 = tuple(new_bond[0])
        br_2, c_2 = tuple(new_bond[-1])

        rwmol = Chem.RWMol(outmol)
        rwmol.AddBond(c_1, c_2, Chem.BondType.SINGLE)

        rwmol.RemoveAtom(br_2)
        rwmol.RemoveAtom(br_1)

        mol = rwmol.GetMol()

        return mol

    def branched_polymer(self, relax_iterations=100, force_field="MMFF", swap_H=True):
        """
        Constructs a complete **branched polymer** by iteratively attaching
        arm molecules to the core molecule, followed by 3D structure optimization.

        This is the main function for generating the final branched polymer.
        It first merges the initial arm to the core. Then, it repeatedly
        attaches additional arms to any remaining placeholder atoms on the
        growing structure. Finally, it performs a geometry optimization
        using the specified force field, with an option to replace the
        placeholder atoms with hydrogen atoms.

        Parameters
        ----------
        relax_iterations : int, optional
            The **maximum number of iterations** to use for the force field
            relaxation process. A higher number generally leads to more
            thorough minimization but increases computation time. Defaults to 100.

        force_field : str, optional
            The **name of the force field** to employ for geometry optimization.
            Currently supported options are "MMFF" (which defaults to MMFF94)
            and "UFF". Defaults to "MMFF".

        swap_H : bool, optional
            A flag indicating whether the placeholder atoms remaining in the
            final branched polymer structure should be **replaced with hydrogen atoms**.
            * If `True` (default), placeholder atoms will be substituted by hydrogens,
              and the molecule will be optimized.
            * If `False`, the placeholder atoms will remain in the structure,
              and the molecule will be optimized without this substitution.

        Return
        -------
        newMol_H : rdkit.Chem.rdchem.Mol
            An **RDKit Mol object** representing the fully constructed and
            geometry-optimized branched polymer. The structure will have
            either hydrogen atoms or the original placeholder atoms at the
            terminal points, depending on the `swap_H` parameter.
        """

        core = self.core
        arm = self.arm
        atom = self.atom

        res = self.merge_arms(core, arm, str(atom))
        res.UpdatePropertyCache(strict=False)
        AllChem.EmbedMolecule(res)

        nm_ph = count_plholder(res, str(atom))

        for _ in range(int(nm_ph)):
            res = self.merge_arms(res, arm, str(atom))

        if swap_H:
            mol = swap_hyd(res, relax_iterations, str(atom), force_field)

        if not swap_H:
            mol = no_swap(res, relax_iterations, force_field)

        return mol


    

