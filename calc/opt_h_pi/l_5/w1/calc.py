import sys
import numpy as np
import pandas as pd
sys.path.append("../../../../src_py/driv_grid")
from driv_grid import solve_driv

sys.path.append("../../../../r1basis")
from r1basis import *
from opt_green import *

channel = '1s->kp'
dipole = 'length'

opt_main(basis_type = 'STO',
         basis_info = 
         [(2, 0.99657  -0.00137j, 'o'),
          (2, 1.00303  -0.28367j, 'o'),
          (2, 0.84629  -0.69526j, 'o'),
          (2, 0.48180  -1.00239j, 'o'),
          (2, 0.14120  -1.06627j, 'o')],
         w0 = 1.0,
         tol = pow(10.0, -10.0),
         target = 'h_pi',
         channel= channel,
         dipole = dipole,
         print_level = 2,
         outfile = "res.out",
         wf_outfile = "wf.csv",
         wf_rs = np.linspace(0, 40.0, 200))

