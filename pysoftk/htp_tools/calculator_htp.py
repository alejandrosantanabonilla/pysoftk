import glob
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor
import os
import fnmatch
from tqdm import tqdm
import time
import shutil

class Htp(object):
    """Class enabling High-throughput quantum chemistry calculations using xtb.

    This class facilitates running xtb calculations (like geometry optimizations)
    on multiple molecular structures organized in subdirectories. It automatically
    scans a specified main directory for subdirectories, checks if each contains
    a single geometry file ('.xyz' by default), and then runs the requested xtb
    calculation method (GFN or GFN-FF) in parallel within each valid subdirectory.

    It intelligently locates the 'xtb' executable, either by searching the system's
    PATH environment variable (common for Conda installations) or by using a
    direct path provided during initialization.

    Attributes
    ----------
    xtb_executable_path : str
        The absolute path to the 'xtb' executable that was found and verified
        during initialization. This path is used for all subsequent calculation calls.
    ending : str
        The default file extension (e.g., '.xyz') used when searching for geometry
        files within subdirectories. This can be customized during initialization
        but defaults to '.xyz' as used by the primary methods.

    Raises
    ------
    FileNotFoundError
        During initialization (__init__) if the specified 'xtb_command' cannot be
        found as an executable in the system PATH, nor resolved as a direct
        valid file path.

    Note
    ----
    - Requires the 'tqdm' library for progress bars (`pip install tqdm`).
    - Requires a working installation of 'xtb' accessible either via the system
      PATH or by providing a direct path to the executable.
    - Assumes a directory structure where each calculation is contained within
      its own subdirectory under a main directory, and each calculation
      subdirectory contains exactly one input geometry file with the expected
      ending (e.g., 'molecule.xyz').

    """

    def __init__(self, xtb_command="xtb", ending=".xyz"):
        """Initializes the Htp class and locates the xtb executable.

        Sets up the instance by determining the default file ending for geometry
        files and, crucially, finding the absolute path to the 'xtb' executable.
        It searches the system's PATH first (using `shutil.which`) and, if not
        found there, checks if the provided `xtb_command` string is a direct
        path to an existing file.

        Parameters
        ----------
        xtb_command : str, optional
            The command name used to invoke the xtb executable (e.g., "xtb")
            which will be searched for in the system PATH, OR the full, direct
            path to the xtb executable file. Defaults to "xtb".
        ending : str, optional
            The default file extension (including the dot) to be associated with
            molecular geometry files for methods that might use it. The primary
            methods currently hardcode '.xyz'. Defaults to ".xyz".

        Raises
        ------
        FileNotFoundError
            If the `xtb_command` cannot be resolved to an actual executable file
            either by searching the PATH or by interpreting it as a direct path.
        """
        # Standardize the file ending format (ensure it starts with '.')
        if not ending.startswith('.'):
            self.ending = "." + ending
        else:
            self.ending = ending

        # Use shutil.which to find the command in PATH
        resolved_path = shutil.which(xtb_command)

        if resolved_path:
            # Found executable in PATH
            self.xtb_executable_path = resolved_path
            print(f"Using xtb executable found at: {self.xtb_executable_path}")
        else:
            # Not found in PATH, check if the input string is a direct file path
            if os.path.isfile(xtb_command):
                # Input is a valid file path
                self.xtb_executable_path = os.path.abspath(xtb_command)
                print(f"Using xtb executable provided at: {self.xtb_executable_path}")
            else:
                # Cannot find the executable by either method
                raise FileNotFoundError(
                    f"Could not find the specified xtb command '{xtb_command}' "
                    f"in PATH, nor does it seem to be a valid file path."
                )

    @staticmethod
    def xtb_gfn(xtb_path, xyz_filename, target_dir, num_cores=None, threshold=None):
        """Executes a single GFN-XTB geometry optimization via subprocess.

        This static method is designed to be called as a target function, often
        in parallel (e.g., via ThreadPoolExecutor). It constructs and runs the
        xtb command line for a GFN-XTB optimization (`--opt`) within the specified
        `target_dir`. Output is redirected to 'output.log' within that directory.

        Parameters
        ----------
        xtb_path : str
            The absolute path to the verified xtb executable.
        xyz_filename : str
            The filename (basename, e.g., "molecule.xyz") of the input geometry
            file. This file is expected to exist within `target_dir`.
        target_dir : str
            The absolute path to the directory where the calculation should run.
            The xtb command will be executed with this directory as its current
            working directory (`cwd`).
        num_cores : int, optional
            The number of processor cores the xtb calculation should use
            (`--parallel` flag). If None, defaults to 1.
        threshold : str, optional
            The convergence threshold level for the geometry optimization
            (`--opt` flag value, e.g., "crude", "normal", "tight"). If None,
            defaults to "crude".

        Returns
        -------
        str
            The absolute path to the output log file created within `target_dir`
            (e.g., "/path/to/target_dir/output.log").

        Side Effects
        ------------
        - Creates or overwrites the file "output.log" in `target_dir`.
        - Runs an external `xtb` process.
        - Prints warning messages to standard output if the `xtb` process
          returns a non-zero exit code or if a Python exception occurs during
          the subprocess execution.
        - Appends error details (Python exception, stderr from xtb) to the
          "output.log" file upon failure.
        """
        
        # Set default values if None
        num_cores = 1 if num_cores is None else int(num_cores)
        threshold = "crude" if threshold is None else str(threshold)

        # Define the output log file path
        out_file = os.path.join(target_dir, "output.log")

        # Construct the command line arguments, quoting paths/filenames
        cmd = (
            f'"{xtb_path}" "{xyz_filename}" '
            f'--parallel {num_cores} '
            f'--opt {threshold} '
        )

        try:
            # Open the output file and run the subprocess
            with open(out_file, 'w') as result_file:
                process = sp.Popen(
                    cmd,
                    stdin=sp.PIPE,       # Standard input pipe (not used by xtb opt)
                    stdout=result_file,  # Redirect standard output to log file
                    stderr=sp.PIPE,      # Capture standard error
                    shell=True,          # Use shell (allows command interpretation, use carefully)
                    cwd=target_dir       # Set working directory for xtb
                )
                # Wait for process to finish and capture stderr
                stderr_output = process.communicate()[1]

                # Check if the process exited with an error code
                if process.returncode != 0:
                    print(f"Warning: xtb calculation in {target_dir} failed.")
                    print(f"Stderr:\n{stderr_output.decode(errors='ignore')}") # Decode stderr
        except Exception as e:
            # Handle Python exceptions during subprocess execution
            print(f"Error running xtb in {target_dir}: {e}")
            # Append Python error and potentially stderr to the log file
            with open(out_file, 'a') as error_file:
                error_file.write(f"\n--- PYTHON ERROR ---\n{e}\n")
                if 'stderr_output' in locals(): # Check if stderr was captured
                    error_file.write(f"\n--- STDERR ---\n{stderr_output.decode(errors='ignore')}\n")
        # Return the path to the log file, regardless of success/failure
        return out_file


    @staticmethod
    def xtb_ff(xtb_path, xyz_filename, target_dir, num_cores=None, threshold=None):
        """Executes a single GFN-FF geometry optimization via subprocess.

        Similar to `xtb_gfn`, but performs a GFN-FF (force field) optimization
        by adding the `--gfnff` flag to the xtb command line. Output is
        redirected to 'output_ff.log' within the `target_dir`.

        Parameters
        ----------
        xtb_path : str
            The absolute path to the verified xtb executable.
        xyz_filename : str
            The filename (basename) of the input geometry file within `target_dir`.
        target_dir : str
            The absolute path to the directory where the calculation should run (cwd).
        num_cores : int, optional
            Number of processor cores for xtb (`--parallel`). Defaults to 1.
        threshold : str, optional
            Convergence threshold for optimization (`--opt`). Defaults to "crude".

        Returns
        -------
        str
            The absolute path to the output log file 'output_ff.log' within
            `target_dir`.

        Side Effects
        ------------
        - Creates or overwrites the file "output_ff.log" in `target_dir`.
        - Runs an external `xtb` process with the `--gfnff` flag.
        - Prints warnings/errors similar to `xtb_gfn`.
        - Appends error details to "output_ff.log" upon failure.
        """
        
        # Set default values if None
        num_cores = 1 if num_cores is None else int(num_cores)
        threshold = "crude" if threshold is None else str(threshold)

        # Define the specific output log file path for GFN-FF
        out_file = os.path.join(target_dir, "output_ff.log")

        # Construct the command line, including the --gfnff flag
        cmd = (
            f'"{xtb_path}" --gfnff "{xyz_filename}" '
            f'--parallel {num_cores} '
            f'--opt {threshold} '
        )

        try:
            # Open the output file and run the subprocess
            with open(out_file, 'w') as result_file:
                process = sp.Popen(
                    cmd,
                    stdin=sp.PIPE,
                    stdout=result_file,
                    stderr=sp.PIPE,
                    shell=True,
                    cwd=target_dir
                )
                # Wait for process and capture stderr
                stderr_output = process.communicate()[1]

                # Check for non-zero exit code
                if process.returncode != 0:
                    print(f"Warning: xtb (gfnff) calculation in {target_dir} failed.")
                    print(f"Stderr:\n{stderr_output.decode(errors='ignore')}")
        except Exception as e:
            # Handle Python exceptions
            print(f"Error running xtb (gfnff) in {target_dir}: {e}")
            # Append error details to the log file
            with open(out_file, 'a') as error_file:
                error_file.write(f"\n--- PYTHON ERROR ---\n{e}\n")
                if 'stderr_output' in locals():
                    error_file.write(f"\n--- STDERR ---\n{stderr_output.decode(errors='ignore')}\n")
        # Return the path to the GFN-FF log file
        return out_file

    def _find_xyz_in_directory(self, directory):
        """Finds all files ending with '.xyz' within a specific directory.

        This is a helper method that performs a non-recursive search for files
        matching the pattern '*.xyz' directly inside the given `directory`.

        Parameters
        ----------
        directory : str
            The path to the directory to search within.

        Returns
        -------
        list[str]
            A list containing the absolute paths of all found '.xyz' files.
            Returns an empty list if no matching files are found.
        """
        
        # Use glob to find files matching the pattern in the specified directory
        xyz_files = glob.glob(os.path.join(directory, '*.xyz'))
        return xyz_files

    def htp_xtb_gfn(self, directory, max_work, num_cores, threshold="crude"):
        """Runs parallel GFN-XTB optimizations on structures in subdirectories.

        This method orchestrates the high-throughput process. It scans the
        specified `directory` for immediate subdirectories. For each subdirectory
        that contains exactly one '.xyz' file, it schedules an `xtb_gfn`
        calculation using a `ThreadPoolExecutor` for parallel execution.
        Progress is displayed using `tqdm`.

        Parameters
        ----------
        directory : str
            The path to the main directory. This directory should contain
            subdirectories, each holding one calculation setup (specifically,
            one '.xyz' file).
        max_work : int
            The maximum number of concurrent `xtb` processes to run in parallel.
            This corresponds to the number of worker threads in the pool.
        num_cores : int
            The number of CPU cores allocated to *each individual* `xtb` process
            (passed via the `--parallel` flag to `xtb_gfn`).
        threshold : str, optional
            The geometry optimization convergence threshold level passed down to
            each `xtb_gfn` call. Defaults to "crude".

        Returns
        -------
        None
            This method primarily executes side effects (running calculations,
            printing output, creating files) and does not return a value.
        """
        
        # --- Scan for Suitable Subdirectories ---
        subdirs_to_process = [] # List to hold paths of valid subdirectories
        xyz_files_map = {}      # Dictionary mapping subdir path to its .xyz file path

        print(f"Scanning directory '{directory}' for subfolders with .xyz files...")
        try:
            # Iterate over items (files and directories) in the main directory
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                # Check if the item is a directory
                if os.path.isdir(item_path):
                    # Find .xyz files specifically within this subdirectory
                    xyz_files = self._find_xyz_in_directory(item_path)

                    # Check if exactly one .xyz file was found
                    if len(xyz_files) == 1:
                        subdirs_to_process.append(item_path)
                        # Store the mapping from directory path to the xyz file path
                        xyz_files_map[item_path] = xyz_files[0]
                    elif len(xyz_files) == 0:
                        # Skip directory if no .xyz file found
                        print(f"  - Skipping '{item}': No .xyz file found.")
                    else:
                        # Skip directory if multiple .xyz files found
                        print(f"  - Skipping '{item}': Found {len(xyz_files)} .xyz files (expected 1).")
        except FileNotFoundError:
            # Handle case where the main directory doesn't exist
            print(f"Error: Main directory '{directory}' not found.")
            return # Stop execution
        except Exception as e:
            # Handle other potential errors during scanning (e.g., permissions)
            print(f"Error scanning directory '{directory}': {e}")
            return # Stop execution

        # Check if any suitable subdirectories were found
        if not subdirs_to_process:
            print("No suitable subdirectories found to process.")
            return # Stop execution
        print(f"\nFound {len(subdirs_to_process)} subdirectories to process.")
        # --- End of Directory Scanning ---


        # --- Parallel Execution using ThreadPoolExecutor ---
        # Create a thread pool with the specified maximum number of workers
        with ThreadPoolExecutor(max_workers=int(max_work)) as executor:
            futures = [] # List to hold Future objects representing submitted tasks
            print(f"Submitting GFN-XTB jobs using up to {max_work} workers...")

            # Iterate through the valid subdirectories found earlier
            for subdir_path in subdirs_to_process:
                # Retrieve the full path to the .xyz file for this subdir
                xyz_full_path = xyz_files_map[subdir_path]
                # Extract just the filename (basename) from the full path
                xyz_filename = os.path.basename(xyz_full_path)

                # Submit the static xtb_gfn method as a task to the executor
                # Pass necessary arguments, including the stored executable path
                futures.append(executor.submit(
                    Htp.xtb_gfn,                      # The function to execute
                    xtb_path=self.xtb_executable_path, # Use stored path from __init__
                    xyz_filename=xyz_filename,        # Pass only the filename
                    target_dir=subdir_path,           # Directory to run in
                    num_cores=int(num_cores),         # Cores per job
                    threshold=str(threshold)          # Optimization threshold
                ))

            # --- Monitor Progress ---
            print("Running calculations...")
            # Use tqdm to iterate over the futures, showing a progress bar
            # as tasks complete (or raise exceptions).
            for future in tqdm(futures, total=len(futures), desc="Systems relaxed (GFN)"):
                try:
                    # future.result() will wait for the task to complete
                    # and return its result (the log file path from xtb_gfn).
                    # It also re-raises any exception that occurred within the task.
                    future.result() # Check for exceptions raised by the task
                except Exception as e:
                    # Catch exceptions from within the xtb_gfn function itself
                    # (though xtb_gfn has its own internal error handling too).
                    print(f"\nAn error occurred processing a task: {e}")
            # --- End of Progress Monitoring ---

        print("All GFN-XTB calculations submitted (or attempted).")


    def htp_xtb_ff(self, directory, max_work, num_cores, threshold="crude"):
        """Runs parallel GFN-FF optimizations on structures in subdirectories.

        This method orchestrates high-throughput GFN-FF calculations. It mirrors
        the functionality of `htp_xtb_gfn` but schedules `xtb_ff` tasks instead,
        utilizing the GFN-FF force field method in xtb.

        Parameters
        ----------
        directory : str
            Path to the main directory containing subdirectories for calculations.
        max_work : int
            Maximum number of concurrent `xtb` processes (worker threads).
        num_cores : int
            Number of CPU cores for *each individual* `xtb` GFN-FF process.
        threshold : str, optional
            Optimization convergence threshold passed to `xtb_ff`. Defaults to "crude".

        Returns
        -------
        None
            Executes calculations and prints progress; does not return a value.
        """
        
        # --- Scan for Suitable Subdirectories (Identical to htp_xtb_gfn) ---
        subdirs_to_process = []
        xyz_files_map = {}
        print(f"Scanning directory '{directory}' for subfolders with .xyz files (for GFN-FF)...")
        try:
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                if os.path.isdir(item_path):
                    xyz_files = self._find_xyz_in_directory(item_path)
                    if len(xyz_files) == 1:
                        subdirs_to_process.append(item_path)
                        xyz_files_map[item_path] = xyz_files[0]
                    elif len(xyz_files) == 0:
                        print(f"  - Skipping '{item}': No .xyz file found.")
                    else:
                        print(f"  - Skipping '{item}': Found {len(xyz_files)} .xyz files (expected 1).")
        except FileNotFoundError:
            print(f"Error: Main directory '{directory}' not found.")
            return
        except Exception as e:
            print(f"Error scanning directory '{directory}': {e}")
            return

        if not subdirs_to_process:
            print("No suitable subdirectories found to process.")
            return
        print(f"\nFound {len(subdirs_to_process)} subdirectories to process for GFN-FF.")
        # --- End of Directory Scanning ---

        # --- Parallel Execution using ThreadPoolExecutor ---
        with ThreadPoolExecutor(max_workers=int(max_work)) as executor:
            futures = []
            print(f"Submitting GFN-FF jobs using up to {max_work} workers...")

            # Iterate through the valid subdirectories
            for subdir_path in subdirs_to_process:
                xyz_full_path = xyz_files_map[subdir_path]
                xyz_filename = os.path.basename(xyz_full_path)

                # Submit the static xtb_ff method as a task
                futures.append(executor.submit(
                    Htp.xtb_ff,                       # Target function for GFN-FF
                    xtb_path=self.xtb_executable_path, # Use stored path
                    xyz_filename=xyz_filename,
                    target_dir=subdir_path,
                    num_cores=int(num_cores),
                    threshold=str(threshold)
                ))

            # --- Monitor Progress (Identical structure to htp_xtb_gfn) ---
            print("Running GFN-FF calculations...")
            for future in tqdm(futures, total=len(futures), desc="Systems relaxed (GFN-FF)"):
                try:
                    future.result() # Wait for completion and check for exceptions
                except Exception as e:
                    print(f"\nAn error occurred processing a GFN-FF task: {e}")
            # --- End of Progress Monitoring ---

        print("All GFN-FF calculations submitted (or attempted).")


