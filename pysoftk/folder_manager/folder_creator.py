import numpy as np
import pathlib
import glob
import os
import sys
import shutil
import uuid
from pathos.pools import ProcessPool

class Fld:
    """
    A utility class for automated folder and file management.

    This class provides methods to create, organize, and manipulate folders
    and files within the current working directory. It can be used as a
    stand-alone application or integrated with other Python projects like PySoftK.
    """

    def __init__(self):
        """
        Initializes the Fld class.

        This constructor currently performs no specific actions but is included
        for future extensibility.
        """
        
        pass

    def fxd_name(self, testname, times):
        """
        Generates an array of sequentially numbered folder names based on a base name.

        This method is useful for creating a series of folders with a consistent
        naming convention, such as 'experiment_1', 'experiment_2', etc.

        Parameters
        ----------
        testname : str
            The base name to be used for the folders (e.g., "my_experiment").
        times : int
            The number of folder names to generate.

        Returns
        -------
        np.ndarray
            An array of strings, where each string is a unique folder name
            (e.g., ["my_experiment_0", "my_experiment_1", ...]).
        """
        
        return np.array([str(testname) + "_" + str(i) for i in range(times)])

    def _make_dir(self, dir_names):
        """
        Creates a single folder in the current working directory.

        This is a private helper method called by other functions within the class
        to perform the actual directory creation.

        Parameters
        ----------
        dir_names : str
            The name of the folder to be created.

        Returns
        -------
        None
            This method creates a folder at the specified path and does not return a value.

        Raises
        ------
        OSError
            If the directory cannot be created (e.g., due to insufficient permissions,
            or if a file with the same name already exists).
        """
        
        dir_cwd = pathlib.Path().absolute()
        os.mkdir("".join((str(dir_cwd), "/", dir_names)))

    def create(self, times=None, fixed_names=None):
        """
        Creates one or more folders in the current working directory.

        You can either specify the number of folders to create (which will be
        assigned unique random names) or provide a pre-defined array of names.

        Parameters
        ----------
        times : int, optional
            The number of folders to create. If provided, unique random names
            will be generated for each folder. Cannot be used with `fixed_names`.
        fixed_names : np.ndarray, optional
            An array of strings, where each string is a desired name for a folder.
            If provided, `times` will be ignored. Cannot be used with `times`.

        Returns
        -------
        None
            This method creates folders and prints a success message to the console.

        Raises
        ------
        ValueError
            If neither `times` nor `fixed_names` is provided, or if both are provided.
        """
        
        if times is None and fixed_names is None:
            raise ValueError("Either 'times' or 'fixed_names' must be provided.")
        if times is not None and fixed_names is not None:
            raise ValueError("Only one of 'times' or 'fixed_names' can be provided.")

        if fixed_names is not None:
            dir_names = fixed_names
        else:
            times = int(0) if times is None else int(times)
            dir_names = np.array([self._unique_name() for i in range(times)])

        list(map(self._make_dir, dir_names))
        print("Successfully created: " + str(len(dir_names)) + " folders")

    def _unique_name(self):
        """
        Generates a globally unique hexadecimal string.

        This private helper method is used to create unique names for folders
        when `create` is called without `fixed_names`.

        Returns
        -------
        str
            A 32-character hexadecimal string representing a unique UUID.
        """
        
        name = uuid.uuid4().hex
        return name

    def seek_files(self, format_extension):
        """
        Searches for files with a specific extension in the current working directory.

        Parameters
        ----------
        format_extension : str
            The file extension to search for (e.g., "txt", "csv", "pdf").
            The leading dot is optional (e.g., both "txt" and ".txt" are valid).

        Returns
        -------
        list[pathlib.Path]
            A list of `pathlib.Path` objects, each representing a file that matches
            the specified extension in the current directory.
        """
        
        query = "".join(("*", ".", str(format_extension).lstrip('.')))
        inp_name = list(pathlib.Path().absolute().glob(query))
        return inp_name

    def _seek_dir(self):
        """
        Identifies and lists all subdirectories within the current working directory.

        This private helper method is used to get a list of existing directories,
        excluding the '__pycache__' directory.

        Returns
        -------
        List[str]
            A sorted list of absolute paths to all subdirectories found in the
            current working directory, excluding '__pycache__'.
        """
        
        directory = pathlib.Path().absolute()
        folder_dir = [f.path for f in os.scandir(directory) if f.is_dir()]
        x = [i for i, item in enumerate(folder_dir) if item.endswith('__pycache__')]
        if x:
            folder_dir.pop(x[0])
        return sorted(folder_dir)

    def copy_dir(self, source, destination):
        """
        Moves a file or directory from a source path to a destination path.

        Note: Despite the name `copy_dir`, this function performs a move operation,
        renaming the source to the destination. If you intend to copy, use `shutil.copy`
        or `shutil.copytree` directly.

        Parameters
        ----------
        source : str
            The path to the file or directory to be moved.
        destination : str
            The new path, including the desired new name, for the moved file or directory.

        Returns
        -------
        None
            This method moves the specified item and does not return a value.

        Raises
        ------
        OSError
            If the move operation fails (e.g., file not found, permission denied).
        """
        
        from pathlib import Path
        Path(source).rename(destination)

    def file_to_dir(self, format_extension, num_cores=None, fixed_names=None):
        """
        Moves files with a specific extension into newly created or existing directories
        in a parallelized manner.

        This method first identifies all files matching the given `format_extension`.
        Then, it either creates new uniquely named folders (if `fixed_names` is not
        provided) or uses the `fixed_names` to create folders. Finally, it moves
        each identified file into a corresponding folder, distributing the work
        across multiple CPU cores for efficiency.

        Parameters
        ----------
        format_extension : str
            The extension of the files to be moved (e.g., "txt", "jpg").
        num_cores : int, optional
            The number of CPU cores to utilize for parallel processing.
            If `None`, the process will default to 1 core (sequential execution).
        fixed_names : np.ndarray, optional
            An array of strings to use as names for the newly created folders.
            If provided, the number of folders created will match the length of
            this array. If `None`, folders will be created with unique random names,
            matching the number of files found.

        Returns
        -------
        None
            This method moves files and prints a message indicating completion.

        Raises
        ------
        NotImplementedError
            This error is currently listed in the docstring but not explicitly
            raised by the code. It might be a placeholder for future validation
            or error handling related to folder creation issues.
        """
        
        import os.path
        from pathos.pools import ProcessPool

        num_cores = int(1) if num_cores is None else int(num_cores)
        files = self.seek_files(format_extension)
        names = [os.path.basename(i) for i in files]

        # Use fixed names if provided, otherwise create folders based on number of files
        if fixed_names is not None:
            self.create(fixed_names=fixed_names)
        else:
            self.create(len(names))

        dirs = self._seek_dir()
        destinations = ["".join((dirs[i], "/", names[i])) for i in range(len(names))]

        with ProcessPool(nodes=int(num_cores)) as pool:
            pool.map(self.copy_dir, files, destinations)

    def move_files_to_folder(self, file_extension, folder_name, num_cores=None):
        """
        Moves all files with a specified extension from the current directory
        into a single, newly created or existing folder. The operation can be
        parallelized across multiple CPU cores.

        Parameters
        ----------
        file_extension : str
            The extension of the files to move (e.g., "txt", "csv", "pdf").
            The leading dot is optional (e.g., "txt" and ".txt" are both valid).
        folder_name : str
            The name of the destination folder. If the folder does not exist,
            it will be created.
        num_cores : int, optional
            The number of CPU cores to use for parallel processing. If `None`,
            the function will automatically detect and use the number of
            available cores.

        Returns
        -------
        None
            This method performs file movements and does not return a value.
        """

        if not file_extension.startswith("."):
            file_extension = "." + file_extension

        # Use os.path.join for platform-independent path construction
        folder_path = os.path.join(pathlib.Path().absolute(), folder_name)

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        files = [f for f in os.listdir('.') if f.endswith(file_extension)]
        destinations = [os.path.join(folder_path, f) for f in files]

        num_cores = os.cpu_count() if num_cores is None else num_cores

        with ProcessPool(nodes=int(num_cores)) as pool:
            pool.map(self.copy_dir, files, destinations)
