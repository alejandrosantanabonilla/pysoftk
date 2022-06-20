from setuptools import setup, find_packages

setup(
    name='pysoftk',
    version='0.1.0',
    author='A. Santana-Bonilla',
    author_email='k2031560@kcl.ac.uk',
    packages= find_packages(exclude=["tests", "tests.*"]),
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    #url='http://pypi.python.org/pypi/TowelStuff/',
    #license='MIT',
    description='PySoftK: Python Soft-Matter Kings College London',
    #long_description=open('README.txt').read(),
    install_requires=[
        "pytest >= 2.7.2",
        "networkx >= 2.7",
        "rdkit-pypi >= 2021.3.1.2",
        "MDAnalysis >= 2.0.0",
        "tqdm >= 4.64.0",
        "pathos >= 0.2.8",
        "pyscf >= 2.0.1",
        "pyscf[geomopt]",
        "pyscf-semiempirical >= 0.1.0",
    ],
)
