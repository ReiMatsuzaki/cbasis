import sys
import numpy as np
import pandas as pd
from itertools import product, combinations

sys.path.append("../../../../r1basis")
from r1basis import *
from opt_green import *

def et_zeta(num, z0, z1):
    r = (z1/z0)**(1.0/(num-1))
    return [z0*r**n for n in range(num)]

## this parameter is used with fitted GTOs in H2plus calculation.
et_zeta = et_zeta(15, 0.01, 30.0)

xos = [0.0001, 0.0002, 0.0003, 0.0005, 0.0007,
       0.001,  0.002,  0.003,  0.005,  0.007]
yos = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3]
#xos = [0.001, 0.009]
#yos = [0.001, 0.03]
zos = [x-1.0j*y for (x,y) in product(xos, yos)]

f = open('res.out', 'w')
sys.stdout = f

for (z0, z1) in combinations(zos, 2):
    basis_info = ([(2, z0, 'o')] +
                  [(2, z1, 'o')] +
                  [(2, z, 'f') for z in et_zeta])

    opt_main(basis_type = 'GTO',
             basis_info = basis_info,
             w0 = 1.0,
             tol = pow(10.0, -5.0),
             target = 'h_pi',
             channel= '1s->kp',
             dipole = 'length',
             print_level = 0)

f.close()    

