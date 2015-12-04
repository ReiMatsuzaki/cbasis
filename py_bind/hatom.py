from linspace import *
from set_l2func import *

import numpy as np

class HAtom:
    def __init__(self, z):
        self.z = z

    def eigenfunc(self, n, l):
        if(n == 1 and l == 0):
            return STO(2.0, 1, 1.0)
        if(n == 2 and l == 0):
            return STO(1.0/np.sqrt(2.0), 1, 0.5) + STO(-1.0/(2.0*np.sqrt(2.0)), 2, 0.5)
        if(n == 2 and l == 1):
            return STO(1.0/(np.sqrt(6.0)*2.0), 2, 0.5)
        else:
            raise("not implemented yet")

    def eigenenergy(self, n):
        return -1.0/(2.0*n*n)

    def length(self, n, l0, l1):
        if(n == 1 and l0 == 0 and l1 == 1):
            return STO(2.0, 2, 1.0)
        else:
            raise("not implemented yet")

    def hop(self, l):
        h0 = -0.5*D2() + (-1.0)*Rm(-1)
        if(l == 0):
            return h0
        else:
            return h0 + l*(l+1)*0.5*Rm(-2)

    def h_minus_ene_op(self, n, l, ene):
        return h0 - ene * Rm(0)
