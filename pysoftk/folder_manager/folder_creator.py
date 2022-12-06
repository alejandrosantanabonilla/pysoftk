import numpy as np
import pathlib 

import glob
import os
import sys

import shutil
import uuid


class Fld:
    """Create folders automatically to be used as stand-alone application or along PySoftK.

    Examples
    --------

    """

    def __init__(self):
        pass
            
    def fxd_name(self, testname):
        """Create an array of fixed names for folders.

        Parameters
        ----------

        testname: str
            Base name to be used as name

        Returns
        -------

        np.ndarray
           An array of names for folders
        """
        return np.array([str(testname)+ "_" +
                         str(i) for i in range(self.times)])


    def _make_dir(self, dir_names):
        """Function to create a folder in the current working directory.

        Parameters
        ----------

        dir_names : str
           Name of the folder that will be created.       

        Returns
        -------

        None
           Creates a folder with a provided name.
        """
        dir_cwd = pathlib.Path().absolute()
        os.mkdir("".join((str(dir_cwd),"/", dir_names)))

        
    def create(self, times=None):
        """Function to create a folder in the current working directory.

        Parameters
        ----------

        times : boolean, optional
           Number of times that a folder will be created       
           
        Returns
        -------

        None : None
           Creates a folder with a provided name.
        """

        times = int(0) if times is None else int(times)
        
        dir_names = np.array([self._unique_name()
                              for i in range(times)])
       
        list(map(self._make_dir,dir_names)) 
        
        print ("Succesfully created: " + str(len(dir_names)) + " folders")
        

    def _unique_name(self):
        """Function to create an unique name
      
        Returns
        -------

        name : str
           Creates a folder with a random name.
        """
        name = uuid.uuid4().hex
        return name


    def seek_files(self, format_extension):
        """Function to seek files in the current working directory.
       
        Parameters
        ----------

        format_extension : str
            Extension used to seek in the current working directory
 
        Returns
        -------

        inp_name : str
            
        """
        query="".join(("*",".",str(format_extension)))
        inp_name = list(pathlib.Path().absolute().glob(query))
        return inp_name


    def _seek_dir(self):
        """Function to seek and list directories in the current working directory.

        Returns
        -------

        folder_dir : List[n str]
           Sorted list of folders inside the current working directory.
        """
        import re
        directory = pathlib.Path().absolute()
        folder_dir = [f.path for f in os.scandir(directory) if f.is_dir()]
        x=[i for i, item in enumerate(folder_dir)
           if item.endswith('__pycache__')]
        if x:
            folder_dir.pop(x[0])      
        return sorted(folder_dir)

    
    def copy_dir(self, source, destination):
        """Function to copy files to directories.

        Parameters
        ----------

        source : str
          Source path where the files is located.

        destination : str
          Source path where the directory is located.
        
        Returns
        -------

        None
           Creates a folder with a provided name.
        """
        from pathlib import Path
        Path(source).rename(destination)

        
    def file_to_dir(self, format_extension, num_cores=None):
        """Function to move files to directories in parallel
  
        Parameters
        ----------

        format_extension : str
            Format extension used to seek.
    
        num_cores : int. optional
            Number of cores to be used.

        Results
        -------

        None:
            Move files to directories.
      
        Raises
        ------
        NotImplementedError
            Folders can not be created.
        """
        import os.path
        
        from pathos.pools import ProcessPool
        
        
        num_cores = int(1) if num_cores is None else int(num_cores)
        
        
        files = self.seek_files(format_extension)
        names=[os.path.basename(i) for i in files]
        
        self.create(len(names))

        dirs = self._seek_dir()
        destinations= ["".join((dirs[i],"/",names[i]))
                       for i in range(len(names))]

        with ProcessPool(nodes=int(num_cores)) as pool:
           pool.map(self.copy_dir, files, destinations)
           pool.close()
           pool.join()
                
