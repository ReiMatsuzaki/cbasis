import sys
import numpy as np
import pandas as pd

sys.path.append("../../../r1basis")
from r1basis import *
from opt_green import *

f = open("res_one.out", "w")
def run():
    h_pi = H_Photoionization(1.0, 1, 0, 1, "length")
    z0s = [0.6-0.6j]
    opt = [True for z in z0s]
    basis = STOs().add(2, z0s).setup()
    res = newton(vgh_green(h_pi, basis, opt), z0s)

    f.write(str(res))

run()    
