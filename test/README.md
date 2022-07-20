# Test-suite

This folder contains the PySoftK unit-test suite. It is written
to be used in combination with the pytest library, which is required and 
installed in the setup.py file. 

# How to run the test-suite?

To run the test-suite, you need to go to the test folder and type

pytest

Once this is done the following result should appear:

# Result of the tests:

============= test session starts ================
platform linux -- Python 3.8.10, pytest-7.1.2, pluggy-1.0.0

plugins: cov-3.0.0

collected 19 items                                                                                                                         

test_chains_func.py .....                                                                                                            [ 26%]

test_folders_func.py ...                                                                                                             [ 42%]

test_format_print.py ....                                                                                                            [ 63%]

test_htp_ff.py .                                                                                                                     [ 68%]

test_htp_gfn.py .                                                                                                                    [ 73%]

test_sm_func.py .....                                                                                                                [100%]

=========== 19 passed in 64.44s (0:01:04) ========
