import os
from pysoftk.tools.utils_ob import *

test_file = os.path.join(os.getcwd(), 'rrwtwe_ba_kk_psk.pdb')
check_bond_order(test_file)
