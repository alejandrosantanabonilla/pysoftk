![PySoftK testing](https://github.com/alejandrosantanabonilla/pysoftk/actions/workflows/pyci.yml/badge.svg)
<!-- Pytest Coverage Comment:Begin -->

<!-- Pytest Coverage Comment:End -->

<p align="center">
    <img width="450" src="pysoftklogo.png" >
</p>


**PySoftK** is a set of Python tools and programs for modelling and simulating polymers with different topologies. The program is still under active 
development and contributions are welcome. A complete introduction into the program can be found in this link [Documentation][1]. 


# Pip installation

To quickly install **PySoftk**, we encourage to do it inside a virtual environment, which can be achieved in the following way:

1. Create a directory named as you want and access it (in this case called work_pol):

```bash 
  [~] mkdir work_pol
  [~] cd work_pol
```

2. Create a virtual environment named pol (this name can be changed, of course) and activate the environment:

```bash 
   [~] python -m venv pol
   [~] source pol/bin/activate
```

3. Get PySoftK from the GitHub repository:

```bash 
  git clone https://github.com/alejandrosantanabonilla/pysoftk.git
```

4. Install PySoftK in this folder using the virtual environment

```bash 
   [~] cd pysoftk
   [~] pip install .
```

**NOTE:** To use the **calculators** option, the code [xtb][2] needs to be installed and linked to the executable using the command:

```bash  
   export PATH="$PATH:$HOME/PATH_TO_FOLDER/xtb-X.X.X/bin"
   export XTBHOME="PATH_TO_XTB_FOLDER/xtb-X.X.X"
   ulimit -s unlimited
```


# Conda installation

These are instructions to 

1. It is recommended to create a new virtual environment in which the [PySoftK][1] dependencies can be installed, for example

```console
   conda create --name pysoftk_env
   conda activate pysoftk_env
```

2. Make a directory to download [PySoftK][1] and change to this directory

```console
   (pysoftk_env) mkdir pysoftk_dir
   (pysoftk_env) cd pysoftk_dir
```

3. Download [PySoftK][1] from Github

```console
   (pysoftk_env) git clone https://github.com/alejandrosantanabonilla/pysoftk.git
```

4. Install [PySoftK][1] and its dependencies as required

```console

   (pysoftk_env) cd pysoftk
   (pysoftk_env) pip install .
```

5. Install [xtb][2] using conda

```console
   conda install -c conda-forge xtb
```

6. You may need to add the [PySoftK][1] installation directory to your .bashrc or .zshrc file

```console 
    PATH=$PATH:/<path to>/pysoftk_dir	
```

That's it, you're ready to use [PySoftK][1] (remember to activate pysoftk_env!)


Replacing **PATH_TO_FOLDER** by the actual path where xtb is stored in your computer.
  
# Testing your installation:
  
5. For testing PySoftK, you need to go to the folder **test** and then type:

```bash 
  pytest
```


[1]: https://alejandrosantanabonilla.github.io/pysoftk/
[2]: https://github.com/grimme-lab/xtb

