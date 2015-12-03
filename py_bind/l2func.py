from l2func_bind import *

class OpAdd:
    def __init__(self, a, b):
        self.left = a
        self.right = b

class OpMult:
    def __init__(self, a, b):
        self.left = a
        self.right = b

cip_dict = {}
cip_dict[(STO, STO)] = cip_ss
cip_dict[(GTO, STO)] = cip_gs
cip_dict[(STO, GTO)] = cip_sg
cip_dict[(GTO, GTO)] = cip_gg

def get_cip(a, b):
    return cip_dict[(type(a), type(b))]

cip_op_dict = {}
cip_op_dict[(STO, OpD2, STO)] = cip_s_d2_s
cip_op_dict[(STO, OpD2, GTO)] = cip_s_d2_g   
cip_op_dict[(GTO, OpD2, STO)] = cip_g_d2_s
cip_op_dict[(GTO, OpD2, GTO)] = cip_g_d2_g
                               
cip_op_dict[(STO, OpD1, STO)] = cip_s_d1_s
cip_op_dict[(STO, OpD1, GTO)] = cip_s_d1_g
cip_op_dict[(GTO, OpD1, STO)] = cip_g_d1_s
cip_op_dict[(GTO, OpD1, GTO)] = cip_g_d1_g
                               
cip_op_dict[(STO, OpRm, STO)] = cip_s_rm_s
cip_op_dict[(STO, OpRm, GTO)] = cip_s_rm_g
cip_op_dict[(GTO, OpRm, STO)] = cip_g_rm_s
cip_op_dict[(GTO, OpRm, GTO)] = cip_g_rm_g

def get_cip_op(a, o, b):
    return cip_op_dict[(type(a), type(o), type(b))]

def cip_a_op_b(a, o, b):
    if(isinstance(o, OpAdd)):
        return cip_a_op_b(a, o.left, b) + cip_a_op_b(a, o.right, b)
    if(isinstance(o, OpMult)):
        return o.left * cip_a_op_b(a, o.right, b)
    return get_cip_op(a, o, b)(a, o, b)

def cip(a, b, c=None):
    if(c == None):
        return get_cip(a, b)(a, b)
    else:
        return cip_a_op_b(a, b, c)



