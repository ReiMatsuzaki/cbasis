import numpy as np
from symmolint import *
def SymGTOs_at_r_ylm(gtos, L, M, irrep, cs, rs):
    (ys, ds) = gtos.at_r_ylm_cpp(L, M, irrep, cs, rs)
    return (np.array(ys), np.array(ds))
    #return np.array(gtos.at_r_ylm_cpp(L, M, irrep, cs, rs))

SymGTOs.at_r_ylm = SymGTOs_at_r_ylm

def calc_v_mat(g1, g2, xyz, q):
    g1p = g1.clone()
    g2p = g2.clone()    
    xyzq = [[xyz[0]], [xyz[1]], [xyz[2]], [q]]
    g1p.set_atoms(xyzq)
    g1p.setup()
    g2p.set_atoms(xyzq)
    g2p.setup()
    mat = calc_mat(g1p, g2p, True)
    return mat["v"]
