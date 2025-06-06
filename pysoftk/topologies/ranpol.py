import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import random
import copy

from pysoftk.linear_polymer import linear_polymer as lp
from pysoftk.linear_polymer import super_monomer as sm

from pysoftk.tools.utils_rdkit import *
from pysoftk.tools.utils_func import *

class Rnp():
    """
    A class for the **computational synthesis of random copolymers** from a
    given set of RDKit monomer molecules.

    This class enables the construction of linear polymer chains where the
    sequence of different monomer types is determined stochastically based
    on user-defined probabilities. It supports both two-component (A-B)
    and three-component (A-B-C) random copolymers, followed by 3D embedding
    and force field optimization.

    Examples
    ---------
    >>> from rdkit import Chem
    >>> from pysoftk.topologies.random import Rnp
    >>> # Define two monomers with placeholder atoms (e.g., [Pt])
    >>> monomer_A = Chem.MolFromSmiles("C([Pt])C")
    >>> monomer_B = Chem.MolFromSmiles("N([Pt])N")
    >>>
    >>> # Initialize the Rnp class for A-B random copolymer
    >>> random_copoly_builder_AB = Rnp(ma=monomer_A, mb=monomer_B, atom="Pt")
    >>>
    >>> # Create a random A-B copolymer of length 10 with 70% A probability
    >>> copolymer_AB = random_copoly_builder_AB.random_ab_copolymer(
    ...     len_polymer=10, pA=0.7, force_field="MMFF"
    ... )
    >>> print(f"Random A-B copolymer atoms: {copolymer_AB.GetNumAtoms()}")
    >>>
    >>> # Define a third monomer (C)
    >>> monomer_C = Chem.MolFromSmiles("O([Pt])O")
    >>> # Create a random A-B-C copolymer of length 15 with pA=0.4, pB=0.3
    >>> copolymer_ABC = random_copoly_builder_AB.random_abc_copolymer(
    ...     mc=monomer_C, len_polymer=15, pA=0.4, pB=0.3, force_field="UFF"
    ... )
    >>> print(f"Random A-B-C copolymer atoms: {copolymer_ABC.GetNumAtoms()}")

    Note
    -----
    The **RDKit package must be installed** in your Python environment for this
    class to function correctly. Ensure that all monomer molecules provided
    contain the same specified placeholder atom for proper bonding.
    This class relies on functionalities from `pysoftk.linear_polymer`
    and `pysoftk.tools`.
    """

    __slots__ = ['ma', 'mb', 'atom']

    def __init__(self, ma, mb, atom):
        """
        Initializes the `Rnp` class with two distinct RDKit monomer molecules
        (monomer A and monomer B) and the common placeholder atom used for linking them.

        These input molecules serve as the fundamental building blocks for the
        random copolymerization process. For three-component copolymers,
        the third monomer (C) is passed directly to the `random_abc_copolymer` method.

        Parameters
        -----------
        ma : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing **monomer A**. This molecule
            must contain a designated placeholder atom.

        mb : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing **monomer B**. This molecule
            must also contain a designated placeholder atom.

        atom : str
            The **atomic symbol** (e.g., "Pt", "Au") of the placeholder atom.
            This specific atom type must be present in all monomer types (`ma`, `mb`, `mc`)
            and will be used to establish bonds during the random polymerization process.
        """

        self.ma = ma
        self.mb = mb
        self.atom = atom

    def random_ab_copolymer(self, len_polymer, pA,
                            relax_iterations=100, force_field="MMFF", swap_H=True):
        """
        Constructs a **random A-B copolymer** of a specified length, where the
        incorporation of monomer A versus monomer B is determined by a user-defined
        probability `pA`. The probability of incorporating monomer B is `1 - pA`.

        This function iteratively builds the polymer chain by randomly selecting
        either monomer A or monomer B at each step based on the given probabilities.
        The final polymer undergoes 3D embedding and force field optimization.

        Parameters
        -----------
        len_polymer : int
            The **total number of monomer units** (A + B) that will form the
            random copolymer chain.

        pA : float
            The **probability (between 0.0 and 1.0)** of incorporating monomer A
            at each step of the polymerization. Consequently, the probability
            of incorporating monomer B is `1 - pA`.

        relax_iterations : int, optional
            The **maximum number of iterations** for the force field relaxation
            algorithm during the final polymer optimization. Defaults to 100.

        force_field : str, optional
            The **name of the force field** to apply during the geometry
            optimization. Accepted values are "MMFF" (which defaults to MMFF94)
            or "UFF". Defaults to "MMFF".

        swap_H : bool, optional
            A flag indicating whether the placeholder atoms (specified by `self.atom`)
            at the termini of the constructed polymer should be **replaced with
            hydrogen atoms**.
            * If `True` (default), placeholder atoms are converted to hydrogens
              before the final optimization.
            * If `False`, the placeholder atoms remain in the structure.

        Return
        -------
        mol : rdkit.Chem.rdchem.Mol
            An **RDKit Mol object** representing the fully constructed and
            geometry-optimized random A-B copolymer.
        """

        ma = self.ma
        mb = self.mb
        atom = self.atom

        rand = np.random.rand()

        if rand < pA:
            m1 = ma
        else:
            m1 = mb

        rand = np.random.rand()

        if rand < pA:
            m2 = ma
        else:
            m2 = mb

        monomer = sm.Sm(m1, m2, str(atom))

        for i in range(len_polymer - 1):
            rand = np.random.rand()

            if rand < pA:
                m3 = ma
            else:
                m3 = mb

            monomer = sm.Sm(monomer.mon_to_poly(), m3, str(atom))

        mol = monomer.mon_to_poly()

        if swap_H:
            newMol_H = swap_hyd(mol, relax_iterations, str(atom), force_field)

        if not swap_H:
            newMol_H = no_swap(mol, relax_iterations, force_field)

        return newMol_H

    def random_abc_copolymer(self, mc, len_polymer, pA,
                             pB, relax_iterations=100, force_field="MMFF", swap_H=True):
        """
        Constructs a **random A-B-C terpolymer** of a specified length, where the
        incorporation of monomer A, B, or C is determined by user-defined
        probabilities `pA` and `pB`. The probability of incorporating monomer C
        is implicitly `1 - pA - pB`.

        This function iteratively builds the polymer chain by randomly selecting
        monomer A, B, or C at each step based on the given probabilities.
        The final polymer undergoes 3D embedding and force field optimization.

        Parameters
        ----------
        mc : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing **monomer C**. This molecule
            must also contain the designated placeholder atom.

        len_polymer : int
            The **total number of monomer units** (A + B + C) that will form the
            random terpolymer chain.

        pA : float
            The **probability (between 0.0 and 1.0)** of incorporating monomer A
            at each step of the polymerization.

        pB : float
            The **probability (between 0.0 and 1.0)** of incorporating monomer B
            at each step of the polymerization. It is expected that `pA + pB <= 1.0`.

        relax_iterations : int, optional
            The **maximum number of iterations** for the force field relaxation
            algorithm during the final polymer optimization. Defaults to 100.

        force_field : str, optional
            The **name of the force field** to apply during the geometry
            optimization. Accepted values are "MMFF" (which defaults to MMFF94)
            or "UFF". Defaults to "MMFF".

        swap_H : bool, optional
            A flag indicating whether the placeholder atoms (specified by `self.atom`)
            at the termini of the constructed polymer should be **replaced with
            hydrogen atoms**.
            * If `True` (default), placeholder atoms are converted to hydrogens
              before the final optimization.
            * If `False`, the placeholder atoms remain in the structure.

        Return
        -------
        mol : rdkit.Chem.rdchem.Mol
            An **RDKit Mol object** representing the fully constructed and
            geometry-optimized random A-B-C terpolymer.
        """

        ma = self.ma
        mb = self.mb
        atom = self.atom

        rand = np.random.rand()
        pAB = pA + pB

        if rand < pA:
            m1 = ma

        elif rand < pAB:
            m1 = mb

        else:
            m1 = mc

        rand = np.random.rand()

        if rand < pA:
            m2 = ma

        elif rand < pAB:
            m2 = mb

        else:
            m2 = mc

        monomer = sm.Sm(m1, m2, str(atom))

        for i in range(len_polymer - 1):
            rand = np.random.rand()

            if rand < pA:
                m3 = ma

            elif rand < pAB:
                m3 = mb

            else:
                m3 = mc

            monomer = sm.Sm(monomer.mon_to_poly(), m3, str(atom))

        mol = monomer.mon_to_poly()

        if swap_H:
            newMol_H = swap_hyd(mol, relax_iterations, str(atom), force_field)
        if not swap_H:
            newMol_H = no_swap(mol, relax_iterations, force_field)

        return newMol_H
