from pysoftk.folder_manager.folder_creator import *

# Test 1:
# Seek files with a given extension like: .smi, .mol, .pbd, etc
#Fld(0).seek_files("smi")

# Test 2: Seek for .mol files and create the needed directories
# using 2 cores.
#Fld().file_to_dir("mol",2)

# Test 3: Create 2 folders with random names using 1 core.
Fld().create(2,1)
