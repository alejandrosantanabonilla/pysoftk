# __init__.py for your chemical utilities package
# Let's assume your package directory is named 'pysoftk.tools'
# So this file would be: pysoftk/tools/__init__.py

# --- Import from utils_func.py ---
# These are general utility and RDKit helper functions.
# Sourced from the utils_func.py file you provided.
from .utils_func import (
    get_file_extension,
    pattern_recon,
    pattern_repl,
    pattern_mol_seq,
    count_plholder,
    atom_neigh,
    tuple_bonds,
    create_pol
)
# If you have other general utility functions in utils_func.py to expose, add them above.

# --- Import from utils_ob.py ---
# These functions are related to Open Babel operations.
from .utils_ob import (
    ff_ob_relaxation,
    rotor_opt,
    global_opt,
    check_bond_order
)

# --- Import from utils_rdkit.py ---
# These functions are *exclusively* for RDKit operations.
from .utils_rdkit import (
    no_swap,
    swap_hyd,
    MMFF_rel,
    UFF_rel,
    plc_holder,
    remove_plcholder,
    etkdgv3_energies
)

# --- Define __all__ for the package's public API ---
# This list specifies what `from pysoftk.tools import *` should import.
# It's good practice to explicitly define this.
__all__ = [
    # Exposed from utils_func.py
    'get_file_extension',
    'pattern_recon',
    'pattern_repl',
    'pattern_mol_seq',
    'count_plholder',
    'atom_neigh',
    'tuple_bonds',
    'create_pol',

    # Exposed from utils_ob.py
    'ff_ob_relaxation',
    'rotor_opt',
    'global_opt',
    'check_bond_order',

    # Exposed from utils_rdkit.py
    'no_swap',
    'swap_hyd',
    'MMFF_rel',
    'UFF_rel',
    'plc_holder',
    'remove_plcholder',
    'etkdgv3_energies',
]

# Optional: A message to confirm package initialization (useful for development)
# print(f"Package 'pysoftk.tools' initialized. Public API: {__all__}")
