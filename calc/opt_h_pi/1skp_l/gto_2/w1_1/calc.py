import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *


z0 = 0.164834921195-0.127041608653j
z1 = 0.905433828185-0.58333509576j

opt_main(basis_type = 'GTO',
         basis_info = [('id', True, 2, z0), ('id', True, 2, z1)],
         w0 = 1.0,
         tol = pow(10.0, -8.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 1,
         wf_outfile = "psi.csv",
         wf_rs = np.linspace(0.0, 40.0, 400),
         outfile = "res.out")

