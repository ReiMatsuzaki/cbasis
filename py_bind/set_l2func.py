from l2func_bind import *
from linspace import cip_dict, cip_op_dict, set_as_func, set_as_op, OpId, op_apply_dict

# ==== Set Func ====
# ---- general ----
map(set_as_func, [GTO, STO])

# ---- repr/str ----
def STO_repr(self):
    return "STO({0},{1},{2})".format(self.c.__repr__(), self.n, self.z.__repr__())

def GTO_repr(self):
    return "GTO({0},{1},{2})".format(self.c.__repr__(), self.n, self.z.__repr__())

def CutSTO_repr(self):
    return "CutSTO({0},{1},{2},{3})".format(self.c.__repr__(), self.n, self.z.__repr__(),
                                            self.r0)

STO.__repr__ = STO_repr
GTO.__repr__ = GTO_repr
CutSTO.__repr__ = CutSTO_repr

# ==== Set Operator ====
# ---- general ----
map(set_as_op, [D1, D2, Rm])

# ---- repr ----
def d1_repr(self):
    return "D1"

def d2_repr(self):
    return "D2"

def rm_repr(self):
    return "Rm({0})".format(self.m)

D1.__repr__ = d1_repr
D2.__repr__ = d2_repr
Rm.__repr__ = rm_repr

# ==== inner product ====
cip_dict[(STO, STO)] = cip_ss
cip_dict[(GTO, STO)] = cip_gs
cip_dict[(STO, GTO)] = cip_sg
cip_dict[(GTO, GTO)] = cip_gg
cip_dict[(CutSTO, CutSTO)] = cip_cut_ss

cip_op_dict[(STO, D2, STO)] = cip_s_d2_s
cip_op_dict[(STO, D2, GTO)] = cip_s_d2_g   
cip_op_dict[(GTO, D2, STO)] = cip_g_d2_s
cip_op_dict[(GTO, D2, GTO)] = cip_g_d2_g
cip_op_dict[(CutSTO, D2, CutSTO)] = cip_cut_s_d2_s
                               
cip_op_dict[(STO, D1, STO)] = cip_s_d1_s
cip_op_dict[(STO, D1, GTO)] = cip_s_d1_g
cip_op_dict[(GTO, D1, STO)] = cip_g_d1_s
cip_op_dict[(GTO, D1, GTO)] = cip_g_d1_g
cip_op_dict[(CutSTO, D1, CutSTO)] = cip_cut_s_d1_s
                               
cip_op_dict[(STO, Rm, STO)] = cip_s_rm_s
cip_op_dict[(STO, Rm, GTO)] = cip_s_rm_g
cip_op_dict[(GTO, Rm, STO)] = cip_g_rm_s
cip_op_dict[(GTO, Rm, GTO)] = cip_g_rm_g
cip_op_dict[(CutSTO, Rm, CutSTO)] = cip_cut_s_rm_s


# ==== op func ====
def d1_sto(dum, sto):
    c = sto.c
    n = sto.n
    z = sto.z
    s1 = STO(c*n,  n-1, z)
    s2 = STO(-z*c, n,   z)
    return s1+s2

def d1_gto(dum, gto):
    c = gto.c
    n = gto.n
    z = gto.z
    s1 = GTO(c*n,      n-1, z)
    s2 = GTO(-2.0*z*c, n+1,   z)
    return s1+s2

op_apply_dict[(D1, STO)] = d1_sto
op_apply_dict[(D1, GTO)] = d1_gto

def rm_sto(rm, s):
    f = STO(s.c, s.n+rm.m, s.z)
    return f

def rm_gto(rm, s):
    f = GTO(s.c, s.n+rm.m, s.z)
    return f
    
op_apply_dict[(Rm, STO)] = rm_sto
op_apply_dict[(Rm, GTO)] = rm_gto

def dr(m):
    if(m==1):
        return D1()
    if(m==2):
        return D2()
    else:
        raise Exception("m=1 or 2")

def rm(m):
    return Rm(m)
