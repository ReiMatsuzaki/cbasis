import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *


z0 = 0.0361960723459-0.0271313771956j
z1 = 0.148847470321-0.145461004611j

opt_main(basis_type = 'GTO',
         basis_info = [('id', True, 2, z0), ('id', True, 2, z1)],
         w0 = 1.0,
         tol = pow(10.0, -8.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 1,
         outfile = "res.out",
         wf_outfile="psi.csv",
         wf_rs = np.linspace(0, 40.0, 400))

