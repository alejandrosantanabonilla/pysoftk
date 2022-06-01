# PySoftK
Python program for modelling and simulating polymers. The program is still under active development and contributions are 
welcome. To install pysoftk, we encourage to do it inside a virtual environment, which can be achieved in the following way:

1.) Create a directory named as you want, in this case, we call it, work_pol

mkdir work_pol

the entes that folder:

cd work_pol

2.) Create a virtual environment named pol (this name can be change, of course)

python -m venv pol
source pol/bin/activate

3.) Get pysfotk from the GitHub repository

git clone https://github.com/alejandrosantanabonilla/pysoftk.git


4.) Install pysoftk in this folder using the virtual environment

pip install -e .

This should work. To test pysoftk (for the moment), you can go to the directory test and type

python ANY_FILE_HERE 

that should work indicating that pysoftk is working. 
