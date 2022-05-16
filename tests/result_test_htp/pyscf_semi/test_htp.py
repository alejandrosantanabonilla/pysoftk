from pysoftk.htp_tools.calculator_htp import *

# High-throughput calculations at the gfn-ff level of theory
#Htp("xyz").htp_xtb_ff("xtb",4,1)

# High-throughput calculation at the gnf-xtb level of theory
#Htp("xyz").htp_xtb_gfn("xtb",4,1)

# High-throughput calculation at the PYSCF-semiempirial level of theory
Htp("xyz").htp_pyscf_semi(4, 1)
