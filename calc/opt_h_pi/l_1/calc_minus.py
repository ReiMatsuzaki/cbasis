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

def calc(basis_type, basis_info, **args):
    """ optimize 
    basis_type : str : STOs or GTOs
    basis_info : [(pn :int, zeta :complex, opt_cmd :str|other)]
    .             represent basis functions and its optimization command
    .   opt_cmd can be choosen as follow:
    .      'o'  :  optimize its orbital exponent
    .      'f'  :  fix its orbital exponent
    target : str : calculation target. ex: "h_pi"
    w0 : energy at start of optimization
    ws : list of w for optimization
    tol : tollerance in Newton method.

    -- only if target = h_pi --
    channel : str : "1s->kp" or "2p->ks" etc
    dipole  : str : "length" or "velocity" etc
    """

    if(args["target"] != "h_pi"):
        raise(Exception("only target='h_pi' is supported"))

    ## ==== settings ====
    opt_list = []
    basis = switch(basis_type, {'STO': STOs(), 'GTO': GTOs()})
    for (pn, zeta, cmd) in basis_info:
        opt_list.append(switch(cmd, {'o': True, 'f': False}))
        basis.add(pn, zeta)
    basis.setup()
    
calc("STO", [(2, 0.7-0.3j, 'o')], target='h_pi', channel='1s->kp', dipole='length')
def run():    
    dw = (w1-w0) / (nw-1)
    ws = [w0 + dw * i for i in range(nw)]
    basis = STOs().add(2, z0s).setup()
    zs = z0s
    num_opt = len([1 for o in opt if o])
    
    for w in ws:
        h_pi = H_Photoionization('1s->kp', "length")
        res = newton(vgh_green(h_pi, w, basis, opt), zs, tol=tol)
        zs = res.x
        f.write("w     = {0}\n".format(w))
        for i in range(num_opt):
            f.write("zeta{0} = {1}\n".format(i, res.x[i]))
            
        f.write("alpha = {0}\n".format(res.val))
        f.write("\n")
        
run()    

