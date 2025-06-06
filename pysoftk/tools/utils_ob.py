from .utils_func import get_file_extension
from openbabel import openbabel as ob
from openbabel import pybel as pb

import numpy as np
import os

def ff_ob_relaxation(mol, FF="MMFF94", relax_iterations=100, ff_thr=1.0e-6):
    """Performs geometry optimization on an Open Babel molecule using a specified force field.

    The molecule's coordinates are updated in place.

    Parameters
    ----------
    mol : pybel.Molecule
        The Open Babel molecule (as a Pybel Molecule object) to be optimized.
        The `OBMol` attribute of this object will be used.
    FF : str, optional
        The force field to be used for optimization.
        Supported options include "MMFF94", "UFF", "GAFF".
        Defaults to "MMFF94".
    relax_iterations : int, optional
        The maximum number of iterations for the conjugate gradients algorithm.
        Defaults to 100.
    ff_thr : float, optional
        The convergence criterion (e.g., energy difference or RMS gradient)
        for the force field optimization. Defaults to 1.0e-6.

    Returns
    -------
    pybel.Molecule
        The input molecule with its geometry optimized. The modification is
        done in-place, and the same molecule object is returned.
    """
    # Find and set up the specified force field
    ff = ob.OBForceField.FindForceField(str(FF))
    ff.Setup(mol.OBMol) # Setup the force field with the molecule's OBMol representation

    # Perform conjugate gradients optimization
    ff.ConjugateGradients(int(relax_iterations), float(ff_thr))

    # Update the molecule's coordinates with the optimized geometry
    ff.GetCoordinates(mol.OBMol)

    return mol

def rotor_opt(mol, FF="MMFF94", rot_steps=125):
    """Performs a rotor search (conformational search) on an Open Babel molecule.

    This function attempts to find low-energy conformations by rotating
    rotatable bonds. The molecule's coordinates are updated in place to
    reflect the best conformation found.

    Parameters
    ----------
    mol : pybel.Molecule
        The Open Babel molecule (as a Pybel Molecule object) for which
        to perform the rotor search.
    FF : str, optional
        The force field to be used during the rotor search.
        Supported options include "MMFF94", "UFF", "GAFF".
        Defaults to "MMFF94".
    rot_steps : int, optional
        The number of steps or iterations for the weighted rotor search.
        This influences the thoroughness of the search. Defaults to 125.

    Returns
    -------
    pybel.Molecule
        The input molecule with its conformation optimized by the rotor search.
        The modification is done in-place.
    """
    # Find and set up the specified force field
    ff = ob.OBForceField.FindForceField(str(FF))
    ff.Setup(mol.OBMol)

    # Perform a fast rotor search (preliminary)
    ff.FastRotorSearch()
    
    # Perform a more thorough weighted rotor search
    # The second argument to WeightedRotorSearch is often related to the number of
    # force field operations per cycle. Here it's half of rot_steps (rounded up).
    ff.WeightedRotorSearch(int(rot_steps), int(np.ceil(rot_steps / 2.0)))
    
    # Update the molecule's coordinates with the best conformation found
    ff.GetCoordinates(mol.OBMol)

    return mol

def global_opt(mol, FF="MMFF94", relax_iterations=150, rot_steps=125, ff_thr=1.0e-6):
    """Performs a global geometry optimization on an Open Babel molecule.

    This typically involves a combination of conformational searching (rotor search)
    and local geometry minimizations (e.g., conjugate gradients) to find a
    low-energy structure. The molecule's coordinates are updated in place.

    Parameters
    ----------
    mol : pybel.Molecule
        The Open Babel molecule (as a Pybel Molecule object) to be optimized.
    FF : str, optional
        The force field to be used. Options include "MMFF94", "UFF", "GAFF".
        Defaults to "MMFF94".
    relax_iterations : int, optional
        The maximum number of iterations for conjugate gradients steps.
        Defaults to 150.
    rot_steps : int, optional
        The number of iterations for the weighted rotor search.
        Defaults to 125.
    ff_thr : float, optional
        The convergence criterion for the final conjugate gradients optimization.
        Defaults to 1.0e-6. A less strict criterion (ff_thr * 100) is used
        for the initial relaxation.

    Returns
    -------
    pybel.Molecule
        The input molecule with its geometry globally optimized.
        The modification is done in-place.
    """
    # Find and set up the specified force field
    ff = ob.OBForceField.FindForceField(str(FF))
    ff.Setup(mol.OBMol)

    # Initial conjugate gradients relaxation (potentially with a looser threshold)
    ff.ConjugateGradients(int(relax_iterations), float(ff_thr * 100))
    
    # Perform weighted rotor search for conformational space exploration
    ff.WeightedRotorSearch(int(rot_steps), int(np.ceil(rot_steps / 2.0)))
    
    # Final conjugate gradients relaxation with a stricter threshold
    ff.ConjugateGradients(int(relax_iterations * 100), float(ff_thr)) # Note: relax_iterations is multiplied by 100 here.

    # Update the molecule's coordinates with the optimized geometry
    ff.GetCoordinates(mol.OBMol)

    return mol


def check_bond_order(file_name):
    """Checks and corrects bond orders for a molecule in a given file.

    This function reads a molecule from a file, removes existing hydrogens,
    attempts to connect disconnected parts (if any), perceives bond orders,
    and then re-adds hydrogens. The modified molecule is written to a new
    file with "_new" appended to its original base name.

    Parameters
    ----------
    file_name : str
        The path to the input molecular file. The file format is inferred
        from its extension and must be supported by Open Babel.

    Returns
    -------
    None
        This function does not return a value but writes a new molecular file
        (e.g., if `file_name` is "molecule.mol", output is "molecule_new.mol").
    """
    # Extract base name and extension from the file path
    base_name, file_extension = os.path.splitext(file_name)
    # Get the format string (e.g., "mol", "sdf") from the extension
    format_ext = file_extension.split(".")[-1]
    
    # Read the first molecule from the input file
    # Note: This assumes a single-molecule file or processes only the first molecule.
    inputfile = pb.readfile(str(format_ext), str(file_name))
    mol = next(inputfile) # Get the first molecule object

    # Perform operations on the OBMol object
    mol.OBMol.DeleteHydrogens()      # Remove all hydrogens
    mol.OBMol.ConnectTheDots()       # Attempt to connect separate fragments if appropriate
    mol.OBMol.PerceiveBondOrders()   # Deduce bond orders based on geometry and atom types
    mol.OBMol.AddHydrogens()         # Add hydrogens according to valence rules

    # Construct the output file name
    new_base_name = "_".join([str(base_name), "new"])
    output_file_name = ".".join([new_base_name, str(format_ext)])
    
    # Write the modified molecule to the new file, overwriting if it exists
    output = pb.Outputfile(str(format_ext),
                           output_file_name,
                           overwrite=True)
    output.write(mol)
    output.close()
