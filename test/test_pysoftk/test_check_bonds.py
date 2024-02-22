from pysoftk.tools.utils_ob import *

from openbabel import openbabel as ob
from openbabel import pybel as pb

import pytest
import os

@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

def test_check_pdbs(rootdir):

    test_file = os.path.join(rootdir, 'data/rrwtwe_ba_kk_psk.pdb')
    check_bond_order(test_file)
  
    test_file2 = os.path.join(rootdir, 'data/cys_gg_rrwtwe_psk.pdb')
    check_bond_order(test_file2)

    new_file1=os.path.join(rootdir,'data/rrwtwe_ba_kk_psk_new.pdb')
    new_file2=os.path.join(rootdir,'data/cys_gg_rrwtwe_psk_new.pdb')
    
    assert os.path.exists(new_file1) == True
    assert os.path.exists(new_file2) == True
    
    os.unlink(new_file1)
    os.unlink(new_file2)
