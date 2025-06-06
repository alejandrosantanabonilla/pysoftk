import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import copy

from pysoftk.topologies.diblock import *
from pysoftk.linear_polymer.linear_polymer import *

from pysoftk.tools.utils_func import *
from pysoftk.tools.utils_rdkit import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

class Rn:
    """
    A class dedicated to the **computational synthesis of circular (ring-shaped) polymers**
    from a given RDKit monomer molecule.

    This class provides a robust method to:
    1. Construct a linear polymer chain of a specified length from a monomer.
    2. Identify the terminal placeholder atoms of this linear chain.
    3. Form a new bond between these terminal atoms to cyclize the polymer,
       thereby creating a ring structure.
    4. Perform subsequent 3D embedding and force field optimization to obtain a stable
       and chemically sensible cyclic molecular geometry.

    Examples
    ---------
    >>> from rdkit import Chem
    >>> from pysoftk.topologies.ring import Rn
    >>> # Define a monomer with placeholder atoms (e.g., [Pt])
    >>> monomer_smiles = "C([Pt])C([Pt])"
    >>> monomer = Chem.MolFromSmiles(monomer_smiles)
    >>>
    >>> # Initialize the Rn class with the monomer and the placeholder atom
    >>> ring_builder = Rn(mol=monomer, atom="Pt")
    >>>
    >>> # Create a ring polymer with 5 monomer units
    >>> final_ring_mol = ring_builder.pol_ring(len_polymer=5, force_field="MMFF", relax_iterations=200)
    >>> print(f"Ring polymer atoms: {final_ring_mol.GetNumAtoms()}")
    >>> print(f"Ring polymer bonds: {final_ring_mol.GetNumBonds()}")

    Note
    -----
    The **RDKit package must be installed** in your Python environment for this
    class to function correctly. Ensure that the placeholder atom specified (e.g., "Pt")
    is present at both ends of the monomer unit to allow for linear chain formation
    and subsequent cyclization. This class relies on functionalities from
    `pysoftk.linear_polymer` and `pysoftk.tools`.
    """

    __slots__ = ['mol', 'atom']

    def __init__(self, mol, atom):
        """
        Initializes the `Rn` class by defining the monomer RDKit molecule
        and the placeholder atom used for linking and cyclization.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the monomer unit from which
            the circular polymer will be constructed. This monomer must have
            at least two placeholder atoms (or one placeholder atom that can
            form two bonds) to allow for linear chain extension and subsequent
            ring closure.

        atom : str
            The **atomic symbol** (e.g., "Pt", "Au") of the placeholder atom.
            This specific atom type must be present in the `mol` and will be
            used to establish bonds during the linear polymerization and
            the final cyclization step.
        """
        
        self.mol = mol
        self.atom = atom

    def pol_ring(self, len_polymer=2, force_field="MMFF",
                 relax_iterations=100, shift=1.25, more_iter=10):
        """
        Constructs a **circular (ring-shaped) polymer** from the specified
        monomer by first creating a linear chain and then closing it into a ring.
        The resulting cyclic polymer undergoes comprehensive 3D embedding and
        force field optimization to achieve a stable conformation.

        This function orchestrates the entire ring polymer synthesis process,
        from linear chain formation to the final optimized cyclic structure.

        Parameters
        ----------
        len_polymer : int, optional
            The **number of monomer units** that will comprise the linear chain
            before cyclization. This directly determines the size of the resulting
            ring polymer. Defaults to 2.

        force_field : str, optional
            The **name of the force field** to apply during the geometry
            optimization steps. Valid options are "MMFF" (which defaults to MMFF94)
            or "UFF". Defaults to "MMFF".

        relax_iterations : int, optional
            The **maximum number of iterations** for the force field relaxation
            algorithm during the initial linear polymer formation and the final
            ring optimization. A higher value generally leads to better energy
            minimization. Defaults to 100.

        shift : float, optional
            The **X-axis translation distance** (in Ångstroms) applied to each
            subsequent monomer copy during the initial linear polymer formation.
            This value influences the initial spacing between monomer units.
            Defaults to 1.25.

        more_iter : int, optional
            A **multiplier for additional optimization iterations** during the
            final force field relaxation of the ring polymer. The `relax_iterations`
            value will be multiplied by this number to perform a more extensive
            final geometry optimization, potentially leading to a more stable
            conformation. Defaults to 10.

        Return
        -------
        pol_ring : openbabel.pybel.Molecule
            An **OpenBabel Pybel Molecule object** representing the fully
            constructed and geometry-optimized circular (ring) polymer.
            This object contains the 3D coordinates and bonding information
            of the final cyclic structure.

        Raises
        ------
        ValueError
            If an invalid `force_field` is provided (i.e., not "MMFF" or "UFF").
        """
        
        mol = self.mol
        atom = self.atom

        if force_field not in ("MMFF", "UFF"):
            raise ValueError(f"Invalid force field: {force_field}")

        patt = Lp(mol, str(atom), int(len_polymer), float(shift)).linear_polymer(relax_iterations=int(relax_iterations),
                                                                                 no_att=False)
        proto_ring = patt.write("pdb")

        ring_rdkit = Chem.MolFromPDBBlock(proto_ring)
        atoms = atom_neigh(ring_rdkit, str(atom))

        rwmol = Chem.RWMol(ring_rdkit)

        rwmol.AddBond(atoms[0][1], atoms[1][1],
                      Chem.BondType.SINGLE)

        rwmol.RemoveAtom(atoms[1][0])
        rwmol.RemoveAtom(atoms[0][0])

        final = rwmol.GetMol()
        Chem.SanitizeMol(final)

        newMol_H = Chem.AddHs(final, addCoords=True)
        AllChem.AssignBondOrdersFromTemplate(newMol_H, newMol_H)

        if force_field == "MMFF":
            apply_force_field(newMol_H, "MMFF", relax_iterations * int(more_iter))

        else:
            apply_force_field(newMol_H, "UFF", relax_iterations * int(more_iter))

        pol_ring = check_proto(newMol_H, str(force_field), int(relax_iterations))

        return pol_ring


def apply_force_field(molecule: Chem.Mol, force_field: str, relax_iterations: int) -> None:
    """
    Applies the specified force field (MMFF or UFF) to an RDKit molecule
    to perform geometry relaxation.

    This helper function is used internally to minimize the energy of a
    molecular conformation.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        The RDKit Mol object whose geometry needs to be optimized.
    force_field : str
        The name of the force field to apply ("MMFF" or "UFF").
    relax_iterations : int
        The maximum number of iterations for the force field minimization algorithm.

    Returns
    -------
    None
        This function modifies the input `molecule` in place by updating
        its conformer with the optimized coordinates.
    """

    if force_field == "MMFF":
        MMFF_rel(molecule, relax_iterations)

    else:
        UFF_rel(molecule, relax_iterations)


def check_proto(mol, force_field="MMFF", relax_iterations=100, rot_steps=1, ff_thr=1.0e-6):
    """
    Performs a final, comprehensive geometry optimization and conformational
    search on a proto-polymer or any RDKit molecule using OpenBabel's
    force field capabilities.

    This function is crucial for refining the 3D structure of the generated
    polymer, ensuring it reaches a stable, low-energy conformation. It
    converts the RDKit molecule to an OpenBabel object, applies the chosen
    force field, and can perform a global optimization to explore multiple
    conformers.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit Mol object representing the proto-polymer or molecule
        to be optimized.

    force_field : str, optional
        The name of the force field to use for relaxation. Accepted values
        are "MMFF" (which maps to MMFF94 in OpenBabel) or "UFF". Defaults to "MMFF".

    relax_iterations : int, optional
        The maximum number of iterations for the force field relaxation
        algorithm. Defaults to 100.

    rot_steps : int, optional
        The number of rotational steps to explore during global optimization
        to find different low-energy conformers. Defaults to 1.

    ff_thr : float, optional
        The force field convergence threshold. The optimization will stop
        when the forces on atoms fall below this value. Defaults to 1.0e-6.

    Return
    -------
    newMol_H : openbabel.pybel.Molecule
        An **OpenBabel Pybel Molecule object** representing the
        geometry-optimized molecule.
    """
    
    last_rdkit = Chem.MolToPDBBlock(mol)

    mol_new = pb.readstring('pdb', last_rdkit)

    if force_field == "MMFF":
        ff = pb._forcefields["mmff94"]
        # Relaxation functions from utils in linear_polymer
        #opt_mol=ff_ob_relaxation(mol_new, relax_iterations=int(relax_iterations), ff_thr=1.0e-6)
        opt_mol = global_opt(mol_new, relax_iterations=int(relax_iterations), rot_steps=int(rot_steps), ff_thr=float(ff_thr))

    else:
        ff = pb._forcefields["uff"]
        # Relaxation functions from utils in linear_polymer
        #opt_mol=ff_ob_relaxation(mol_new, relax_iterations=int(relax_iterations), ff_thr=1.0e-6)
        opt_mol = global_opt(mol_new, relax_iterations=int(relax_iterations), rot_steps=int(rot_steps), ff_thr=float(ff_thr))

    return opt_mol
