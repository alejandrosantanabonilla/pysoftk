from setuptools import setup, find_packages

setup(
<<<<<<< HEAD
    name='pysoftk', 
    version='1.0.0', 
    author='A. Santana-Bonilla, R. Lopez-Castro, R. Ziolek, C. Lorenz', 
    author_email='k2031560@kcl.ac.uk', 
=======
    name='pysoftk',
    version='1.0.0',
    author='A. Santana-Bonilla, R. Lopez-Castro, R. Ziolek, C. Lorenz',
    author_email='k2031560@kcl.ac.uk',
>>>>>>> 88f2163339a8017b72d141956fc54ab10ba533f5
    packages= find_packages(exclude=["test", "test.*", "test2", "test.*"]),
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    #url='http://pypi.python.org/pypi/TowelStuff/',
    #license='MIT',
    description='PySoftK: Python Soft-Matter Kings College London',
    #long_description=open('README.txt').read(),
    install_requires=[
        'pytest >= 2.7.2',
<<<<<<< HEAD
        'Cython',
=======
	'Cython',
>>>>>>> 88f2163339a8017b72d141956fc54ab10ba533f5
        'networkx >= 2.7; python_version>="3.8"',
        'networkx < 2.7; python_version<"3.8"',
        'rdkit >= 2021.3.1.2',
        'wheel',
        'tqdm >= 4.64.0',
        'pathos >= 0.2.8',
        'pyscf-semiempirical @ git+https://github.com/pyscf/semiempirical',
        'pyberny',
        'openbabel-wheel',
<<<<<<< HEAD
        'numba',
=======
        'numba < 0.59.0',
>>>>>>> 88f2163339a8017b72d141956fc54ab10ba533f5
	'umap-learn',
	'hdbscan',
	'pandas',
	'pyarrow',
	'MDAnalysis',
    ],
)
