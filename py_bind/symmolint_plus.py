import numpy as np
from symmolint import *
def SymGTOs_at_r_ylm(gtos, L, M, irrep, cs, rs):
    return np.array(gtos.at_r_ylm_cpp(L, M, irrep, cs, rs))

SymGTOs.at_r_ylm = SymGTOs_at_r_ylm
