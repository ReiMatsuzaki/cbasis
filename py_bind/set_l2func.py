from l2func_bind import *
from linspace import cip_dict, cip_op_dict, set_as_func, set_as_op, OpId

# ==== Set Func ====
# ---- general ----
map(set_as_func, [GTO, STO])

# ---- repr/str ----
def STO_repr(self):
    return "STO({0},{1},{2})".format(self.c, self.n, self.z)

def GTO_repr(self):
    return "GTO({0},{1},{2})".format(self.c, self.n, self.z)

STO.__repr__ = STO_repr
GTO.__repr__ = GTO_repr

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

cip_op_dict[(STO, D2, STO)] = cip_s_d2_s
cip_op_dict[(STO, D2, GTO)] = cip_s_d2_g   
cip_op_dict[(GTO, D2, STO)] = cip_g_d2_s
cip_op_dict[(GTO, D2, GTO)] = cip_g_d2_g
                               
cip_op_dict[(STO, D1, STO)] = cip_s_d1_s
cip_op_dict[(STO, D1, GTO)] = cip_s_d1_g
cip_op_dict[(GTO, D1, STO)] = cip_g_d1_s
cip_op_dict[(GTO, D1, GTO)] = cip_g_d1_g
                               
cip_op_dict[(STO, Rm, STO)] = cip_s_rm_s
cip_op_dict[(STO, Rm, GTO)] = cip_s_rm_g
cip_op_dict[(GTO, Rm, STO)] = cip_g_rm_s
cip_op_dict[(GTO, Rm, GTO)] = cip_g_rm_g


