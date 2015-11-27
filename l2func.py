import sympy as sp
import numpy as np
from l2func_bind import *

def sto(c, n, z):
    res = STOs()
    res.add_one(1.0, STO(c, n, z))
    return res

def h_like_atom(s, z=1.0):
    if(s == "1s"):
        return HLikeAtom(1, z, 0)
    elif(s == "2s"):
        return HLikeAtom(2, z, 0)
    elif(s == "2p"):
        return HLikeAtom(2, z, 1)
    elif(s == "3s"):
        return HLikeAtom(3, z, 0)
    elif(s == "3p"):
        return HLikeAtom(3, z, 1)
    else:
        return None

def scalar_prod(scalar, other):
    res = other.clone()
    res.scalar_prod(scalar)
    return res


# add print method for STO and GTO
def exp_basis_show(self, name):
    if self.size() == 0:
        return ''
    def i_th(i):
        c = self.coef_i(i)
        u = self.prim_i(i)
        return '{0}{1}[{2},{3},{4}] + '.format(np.real_if_close(c), 
                                               name, 
                                               np.real_if_close(u.c), 
                                               u.n, 
                                               np.real_if_close(u.z))
    return reduce(lambda x,y:x+y, [ i_th(i) for i in range(self.size()) ])[0:-3]

def sto_show(self):
    return exp_basis_show(self, 'STO')

def gto_show(self):
    return exp_basis_show(self, 'GTO')


def exp_basis_repr(self, name):
    cum_str = ''
    cum_str += name + '\n'
    for i in range(self.size()):
        ui = self.prim_i(i)
        cum_str += str(i) + '\n'
        cum_str += 'd_i: ' + str(self.coef_i(i)) + '\n'
        cum_str += 'c_i: ' + str(ui.c) + '\n'
        cum_str += 'n_i: ' + str(ui.n) + '\n'
        cum_str += 'z_i: ' + str(ui.z) + '\n'
    return cum_str

def sto_repr(self):
    return exp_basis_repr(self, 'STO')

def gto_repr(self):
    return exp_basis_repr(self, 'GTO')


STOs.__str__ = sto_repr
GTOs.__str__ = gto_repr
STOs.__repr__ = sto_show
GTOs.__repr__ = gto_show

# __call__ for operator
Op_sto.__call__ = Op_sto.operate
Op_gto.__call__ = Op_gto.operate

# return function its self as sympy symbol 
def sympy_symbol(m):

    def __func__(self, x):
        cum = 0
        for i in range(self.size()):
            prim = self.prim_i(i)
            cum += self.coef_i(i) * prim.c * (x ** prim.n) * sp.exp(-prim.z * (x ** m))
        return cum
    return __func__

STOs.symbol = sympy_symbol(1)
GTOs.symbol = sympy_symbol(2)

def make_stos(c_sto_list):
    cum = STOs()
    for (c, sto) in c_sto_list:
        cum.add(scalar_prod(c, sto))
    return cum

def make_gtos(c_gto_list):
    cum = GTOs()
    for (c, u) in c_gto_list:
        cum.add(scalar_prod(c, u))
    return cum

def make_op_s(c_op_list):
    cum = Op_sto()
    for (c, o) in c_op_list:
        o.scalar_prod(c)
        cum.add(o)
    return cum

def make_op_g(c_op_list):
    cum = Op_gto()
    for (c, o) in c_op_list:
        o.scalar_prod(c)
        cum.add(o)
    return cum




"""
class STO:
    def sym_ip_with(other):
        return other.sym_ip_with_sto(self)
    def sym_ip_with_gto(gto):
        return sym_ip_sg(self, gto)
    def sym_ip_with_sto(sto):
        return sym_ip_ss(self, sto)

class GTO:
    def sym_ip_with(other):
        return other.sym_ip_with_sto(self)
    def sym_ip_with_gto(gto):
        return sym_ip_gg(self, gto)
    def sym_ip_with_sto(sto):
        return sym_ip_gs(self, sto)    

def sym_ip(a, b):
    return a.sym_ip_with(b)
"""



