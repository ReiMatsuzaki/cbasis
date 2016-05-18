import numpy as np
from r1gtoint_bind import *
def R1GTOs_at_r(gtos, cs, rs):
    ys = gtos.at_r_cpp(cs, rs)
    return np.array(ys)

R1GTOs.at_r = R1GTOs_at_r


def R1GTOs_deriv_at_r(gtos, cs, rs):
    ys = gtos.deriv_at_r_cpp(cs, rs)
    return np.array(ys)

R1GTOs.deriv_at_r = R1GTOs_deriv_at_r


def R1GTOs_deriv_2_at_r(gtos, cs, rs):
    ys = gtos.deriv_2_at_r_cpp(cs, rs)
    return np.array(ys)

R1GTOs.deriv_2_at_r = R1GTOs_deriv_2_at_r


def R1STOs_at_r(stos, rs):
    ys = stos.at_r_cpp(rs)
    return np.array(ys)

R1STOs.at_r = R1STOs_at_r
