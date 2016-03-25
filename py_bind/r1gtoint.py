import numpy as np
from r1gtoint_bind import *
def R1GTOs_at_r(gtos, cs, rs):
    ys = gtos.at_r_cpp(cs, rs)
    return np.array(ys)

R1GTOs.at_r = R1GTOs_at_r
