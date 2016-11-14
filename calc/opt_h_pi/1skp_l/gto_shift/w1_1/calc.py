import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## see old calculation in ~/src/git/opt_cbf/py_bind/calc/opt_gto/ret_shift/
# ==== ET-basis ====
def et_ox_from_max(r0, r1, n, num):
    """
    u(r) = r**n * exp(-zr**2)
    u'(r)= [nr**(n-1) -2zr**(n+1)] * exp(-zr**2)
    <=> n = 2zrr
    <=> z = n/(2* r_max**2)
    """
    z0 = n / (2.0 * r1*r1)
    z1 = n / (2.0 * r0*r0)
    ratio = (z1/z0)**(1.0/(num-1))
    return [z0*ratio**n for n in range(num)]

"""
zs = et_ox_from_max(0.01, 20.0, 2, 20)
print zs[0], zs[-1], len(zs), zs[1]/zs[0]
zs = et_ox_from_max(0.1, 10.0, 2, 20)
print zs[0], zs[-1], len(zs), zs[1]/zs[0]
zs = et_ox_from_max(1.0, 20.0, 2, 10)
print zs[0], zs[-1], len(zs), zs[1]/zs[0]

0.0025 10000.0 20 2.22575235544
0.01 100.0 20 1.62377673919
0.0025 1.0 10 1.94588771758
"""


#z0s = [30.0/(1.7716**n) for n in range(15)]
#z0s = [0.0025*2.225**n for n in range(20)]
#z0s = et_ox_from_max(0.01, 20.0, 2, 20)
#y0 = 0.0137-0.338j
#y0 = 0.0-0.000338j
#y0 = 0.000939964489964-0.0340181892741j
z0s = et_ox_from_max(1.0, 10.0, 2, 10)
#y0 = 0.0137-0.0338j
y0 = 0.0037-0.039j

print z0s[0], z0s[-1], z0s[1]/z0s[0], len(z0s)

basis_info = [('shift', True, 2, z0s, y0)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,

         target = 'h_pi',
         w0 = 1.0,
         channel= '1s->kp',
         dipole = 'length',
         
         tol = pow(10.0, -5.0),
         maxit = 50,
         conv = 'grad',
         
         print_level = 5,
         outfile = "res.out",
         wf_outfile="psi.csv",
         wf_rs = np.linspace(0, 40.0, 400))

