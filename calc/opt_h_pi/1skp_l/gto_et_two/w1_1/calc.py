import sys
import numpy as np
import pandas as pd

## from main_opt_gto.org/Results/FOWF/N-Opt-GTO + M-ET-GTO/2opt-20et

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

basis_info = [('id', True, 2, 0.0016675881-0.058558696j),
              ('id', True, 2, 0.1004093981-0.075346307j),
              ("geo", False, 2, 20, 2.5**(-12), 2.5)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,
         w0 = 1.0,
         tol = pow(10.0, -8.0),
         maxit = 50,
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 2,
         conv = 'dx',
         outfile = "res.out",
         wf_outfile = "psi.csv",
         wf_rs = np.linspace(0.0, 40.0, 200))
