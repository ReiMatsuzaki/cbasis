import sys
import numpy as np
import pandas as pd

sys.path.append("../../../r1basis")
from r1basis import *
from opt_green import *

f = open("res_minus.out", "w")
w0 = 1.0
w1 = 0.5
nw = 11
z0s = [0.6-0.6j]
opt = [True for z in z0s]
tol = pow(10.0, -10.0)

def run():
    dw = (w1-w0) / (nw-1)
    ws = [w0 + dw * i for i in range(nw)]
    basis = STOs().add(2, z0s).setup()
    zs = z0s
    num_opt = len([1 for o in opt if o])
    
    for w in ws:
        h_pi = H_Photoionization(w, 1, 0, 1, "length")
        res = newton(vgh_green(h_pi, basis, opt), zs, tol=tol)
        zs = res.x
        f.write("w     = {0}\n".format(w))
        for i in range(num_opt):
            f.write("zeta{0} = {1}\n".format(i, res.x[i]))
            
        f.write("alpha = {0}\n".format(res.val))
        f.write("\n")
        
run()    

