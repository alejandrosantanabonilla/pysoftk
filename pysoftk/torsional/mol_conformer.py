from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as molDG
from rdkit.Chem import TorsionFingerprints
from pysoftk.tools.utils_rdkit import *
import numpy as np
import os
import sys
from openbabel import pybel, openbabel as ob

class ConformerGenerator:
    """
    A class to generate and save molecular conformers using Open Babel.
    Includes methods for genetic algorithm (GA) based conformer generation
    (modified from the original 'generate_conformers'),
    diversity-based conformer generation (confab), and saving conformers
    to separate files.
    """
    def __init__(self, num_conformers: int = 10, forcefield: str = "mmff94", make3D_steps: int = 50, convergence: int = 5):
        """
        Initializes the ConformerGenerator with default parameters for conformer generation.

        Args:
            num_conformers: The default target number of conformers to generate.
            forcefield: The default force field to use for initial 3D coordinate generation.
            make3D_steps: The default number of steps to use for initial 3D coordinate generation.
            convergence: The default parameter controlling the convergence criteria of the GA-like search.
        """
        self.default_num_conformers = num_conformers
        self.default_forcefield = forcefield
        self.default_make3D_steps = make3D_steps
        self.default_convergence = convergence

    def ga_generate_conformers(self, input_file: str = None,
                                   smiles: str = None,
                                   num_conformers: int = None,
                                   forcefield: str = None,
                                   make3D_steps: int = None,
                                   mutability: int = 5,
                                   convergence: int = None,
                                   num_children: int = None) -> pybel.Molecule:
        """
        Generates conformers using Open Babel's weighted rotor search (similar to a
        genetic algorithm approach). The molecule can be provided as a file
        or a SMILES string.

        Args:
            input_file: Path to a molecular file (e.g., MOL, SDF, XYZ).
                        If provided, this will be used as the input.
            smiles: A SMILES string representing the molecule.
                    This will be used if input_file is not provided.
            num_conformers: The target number of conformers to generate (default: initialized value).
            forcefield: The force field to use for initial 3D coordinate generation
                        (default: initialized value).
            make3D_steps: The number of steps to use for initial 3D coordinate
                          generation (default: initialized value).
            mutability: Parameter controlling the rotor key search randomization
                        (default: 5). Higher values explore more diverse torsions.
            convergence: Parameter controlling the convergence criteria of the
                         search (default: initialized value). Higher values mean stricter convergence.
            num_children: Number of children conformers generated in each step
                          (default: num_conformers * 2).

        Returns:
            A *new* pybel.Molecule object containing the generated conformers.
            The original molecule object is not modified. Returns None if
            the input is invalid or if conformer generation fails.

        Raises:
            ValueError: If num_conformers is less than 1.
            ValueError: If neither input_file nor smiles is provided.
            FileNotFoundError: If input_file is provided but does not exist.
        """
        num_conformers = num_conformers if num_conformers is not None else self.default_num_conformers
        forcefield = forcefield if forcefield is not None else self.default_forcefield
        make3D_steps = make3D_steps if make3D_steps is not None else self.default_make3D_steps
        convergence = convergence if convergence is not None else self.default_convergence

        if num_conformers < 1:
            raise ValueError("Number of conformers (num_conformers) must be at least 1.")

        if input_file:
            if not os.path.exists(input_file):
                raise FileNotFoundError(f"Input file not found: {input_file}")
            try:
                molecule = next(pybel.readfile(os.path.splitext(input_file)[1][1:], input_file))
                base_name, _ = os.path.splitext(os.path.basename(input_file))
                molecule.title = base_name.capitalize()
                print(f"Loaded molecule from: {input_file} ({molecule.title})")
            except Exception as e:
                print(f"Error reading input file '{input_file}': {e}", file=sys.stderr)
                return None
        elif smiles:
            try:
                molecule = pybel.readstring("smi", smiles)
                molecule.title = "Untitled"
                print(f"Loaded molecule from SMILES: {smiles}")
            except Exception as e:
                print(f"Error reading SMILES string '{smiles}': {e}", file=sys.stderr)
                return None
        else:
            raise ValueError("Either 'input_file' or 'smiles' must be provided.")

        # 2. Generate initial 3D coordinates if not present
        if not molecule.OBMol.Has3D():
            print(f"Generating initial 3D coordinates for {molecule.title} using forcefield '{forcefield}'...")
            molecule.make3D(forcefield=forcefield, steps=make3D_steps)
            if not molecule.OBMol.Has3D():
                print("Failed to generate 3D coordinates.", file=sys.stderr)
                return None
            print("Initial 3D coordinates generated.")
        else:
            print(f"Molecule '{molecule.title}' already has 3D coordinates.")

        # Work on a clone to avoid modifying the original molecule
        mol_clone = pybel.Molecule(molecule.OBMol)

        # Determine the number of children if not specified
        children = num_children if num_children is not None else num_conformers * 2

        # --- Conformer Search Setup ---
        cs = ob.OBConformerSearch()
        setup_successful = cs.Setup(mol_clone.OBMol,
                                     num_conformers,   # target number of conformers
                                     children,         # numChildren in docs
                                     mutability,       # mutability in docs
                                     convergence)      # convergence in docs

        if not setup_successful:
            print(f"Error: OBConformerSearch setup failed for molecule: {molecule.title}", file=sys.stderr)
            return None # Or raise an exception

        # --- Perform the Search ---
        print(f"Starting conformer search for {molecule.title}...")
        cs.Search()
        print("Conformer search finished.")

        # --- Retrieve Conformers ---
        cs.GetConformers(mol_clone.OBMol)

        # Get the actual number found (might be less than requested)
        actual_conformers = mol_clone.OBMol.NumConformers()
        print(f"Found {actual_conformers} conformers.")

        return mol_clone

    def confab_generate_conformers(self, molecule_input: str,
                                       output_dir: str,
                                       output_format: str = 'xyz',
                                       base_filename: str = None,
                                       forcefield: str = None,
                                       nconfs: int = 10000,
                                       rmsd: float = 0.5,
                                       energy_gap: float = 50.0,
                                       initial_opt_steps: int = 50,
                                       verbose: bool = False) -> int:
        """
        Generates diverse conformers using a method similar to Confab.

        Args:
            molecule_input: A SMILES string representing the molecule.
            output_dir: The directory to save the generated conformers.
            output_format: The file format for saving (default: 'xyz').
            base_filename: The base name for the output files.
            forcefield: The force field to use for initial 3D generation (default: initialized value).
            nconfs: The maximum number of diverse conformers to generate (default: 10000).
            rmsd: The RMSD threshold for considering conformers as distinct (default: 0.5).
            energy_gap: The energy window (in kJ/mol) to consider for diverse conformers (default: 50.0).
            initial_opt_steps: The number of optimization steps for the initial 3D structure (default: 50).
            verbose: A boolean to control the verbosity (default: False).

        Returns:
            The number of conformers successfully generated and saved.
        """
        forcefield = forcefield if forcefield is not None else self.default_forcefield

        mol = None
        try:
            mol = pybel.readstring("smi", molecule_input)
            mol.addh()
            mol.make3D(forcefield=forcefield, steps=initial_opt_steps)
            ff = pybel._forcefields[forcefield]
            if not ff.Setup(mol.OBMol):
                print(f"Warning: Force field '{forcefield}' setup failed.")
                return 0
            ff.DiverseConfGen(rmsd, nconfs, energy_gap, verbose)
            ff.GetConformers(mol.OBMol)
            num_conformers = mol.OBMol.NumConformers()
            os.makedirs(output_dir, exist_ok=True)
            base = base_filename if base_filename else "molecule"

            # Use the save_conformers_to_separate_files method
            self.save_conformers_to_separate_files(mol, base, output_dir, output_format)

            print(f"Successfully generated and saved {num_conformers} conformers.")
            return num_conformers
        except ValueError as e:
            print(f"Error: {e}")
            return 0
        except KeyError as e:
            print(f"Error: {e}")
            return 0
        except RuntimeError as e:
            print(f"Error: {e}")
            return 0
        except OSError as e:
            print(f"Error: {e}")
            return 0
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return 0

    def save_conformers_to_separate_files(self, molecule_with_conformers: pybel.Molecule,
                                               base_name: str,
                                               output_dir: str = None,
                                               file_format: str = "xyz"):
        """
        Splits a pybel.Molecule object containing multiple conformers into
        separate files, one for each conformer.

        Args:
            molecule_with_conformers: A pybel.Molecule object that contains
                                      multiple conformers (e.g., the output
                                      from generate_conformers()).
            base_name: The base name to use for the output filenames
                       (e.g., "propanol").
            output_dir: The directory where the individual conformer files
                        will be saved. If None, a directory named
                        "{base_name}_conformers" will be created in the
                        current working directory (default: None).
            file_format: The file format to use for saving the conformers
                         (default: "xyz"). Common options include "mol", "pdb", "xyz".
        """
        if not molecule_with_conformers or molecule_with_conformers.OBMol.NumConformers() == 0:
            print("No conformers found in the input molecule, nothing to save.")
            return

        num_found = molecule_with_conformers.OBMol.NumConformers()
        print(f"\nSaving {num_found} conformers for {molecule_with_conformers.title} to separate files.")

        # Determine the output directory
        if output_dir is None:
            output_dir = f"{base_name}_conformers"
        os.makedirs(output_dir, exist_ok=True)
        print(f"Saving individual conformers into directory: '{output_dir}/'")

        # Determine padding for filenames
        padding = len(str(num_found))

        # Loop through each found conformer index
        for i in range(num_found):
            # Set the molecule's coordinates to the i-th conformer
            molecule_with_conformers.OBMol.SetConformer(i)

            # Create a *new* molecule object containing ONLY the current conformer
            single_conformer_mol = pybel.Molecule(molecule_with_conformers.OBMol)

            # Optional: Set a unique title for the individual conformer file
            single_conformer_mol.title = f"{molecule_with_conformers.title}_conformer_{i+1}"

            # Define the output filename
            output_filename = os.path.join(
                output_dir,
                f"{base_name}_conf_{str(i+1).zfill(padding)}.{file_format}"
            )

            # Save the single-conformer molecule to its own file
            try:
                single_conformer_mol.write(file_format, output_filename, overwrite=True)
                print(f" Saved: {os.path.basename(output_filename)}", end='\r')
            except Exception as e:
                print(f"\nError saving conformer {i+1} to file '{output_filename}': {e}", file=sys.stderr)

        print("\nFinished saving all conformers to separate files.")

class Mcon(object):
    """
    A class designed to **generate and manage molecular conformers** for a
    given RDKit molecule.

    This class provides functionality to compute a specified number of 3D
    conformers (different spatial arrangements) for a molecule, leveraging
    RDKit's conformer generation algorithms. The generated conformers,
    along with their calculated energies, can then be saved to a common
    chemical file format (SDF).

    Example
    --------
    >>> from rdkit import Chem
    >>> from pysoftk.conformer.conformer import Mcon
    >>>
    >>> # Create an RDKit molecule from a SMILES string
    >>> mol = Chem.MolFromSmiles("CCCC")
    >>> mol = Chem.AddHs(mol) # Add hydrogens for more realistic conformers
    >>>
    >>> # Initialize the Mcon class to compute 10 conformers
    >>> conformer_generator = Mcon(mol=mol, num_conf=10)
    >>>
    >>> # Generate and save the conformers to an SDF file
    >>> conformer_generator.conformer(output_name="butane_conformers")
    Number of conformers created: 10

    Note
    -----
    This class **requires the RDKit package to be installed** in your Python
    environment. The quality and diversity of the generated conformers
    depend on the underlying RDKit algorithms and the complexity of the
    molecule.
    """

    __slots__ = ['mol', 'num_conf']

    def __init__(self, mol, num_conf):
        """
        Initializes the `Mcon` class with the molecule for which conformers
        are to be generated and the desired number of conformers.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The **RDKit Mol object** representing the molecule for which
            conformers will be calculated. It's often beneficial to add
            hydrogens to the molecule (`Chem.AddHs(mol)`) before passing
            it to this class for more accurate conformer generation.

        num_conf : int
            The **total number of conformers** requested to be computed
            for the provided molecule. The actual number of conformers
            generated might be less if the algorithm cannot find the
            requested amount within its internal limits.
        """
        self.mol = mol
        self.num_conf = num_conf

    def conformer(self, output_name):
        """
        Generates molecular conformers for the initialized molecule using
        RDKit's ETKDGv3 algorithm and saves them to an SDF (Structure-Data File)
        file.

        This method orchestrates the conformer generation process, including
        the calculation of energies for each conformer. It then writes
        all generated conformers into a single SDF file, making them
        accessible for further analysis or visualization. By default, it returns
        all generated conformers, including the highest-energy ones, allowing
        for a comprehensive conformational analysis.

        Parameters
        ----------
        output_name : str
            The **base name** for the output SDF file. The file will be saved
            as `<output_name>.sdf` in the current working directory.

        Return
        -------
        list of rdkit.Chem.Mol objects
            This function doesn't explicitly return a list of molecules directly.
            Instead, it **writes the conformers to an SDF file**.
            The SDF file contains all generated conformers as separate entries,
            each representing an `rdkit.Chem.Mol` object with 3D coordinates.
            A message indicating the number of created conformers is printed
            to the console.
        """
        
        mol = self.mol
        num_conf = self.num_conf

        mol, cids, energies = etkdgv3_energies(mol, num_conf)

        full_name = [output_name, ".sdf"]
        w = AllChem.SDWriter("".join(full_name))

        for cid in cids:
            w.write(mol, confId=cid)
        w.close()

        print("Number of conformers created: {} ".format(len(cids)))
