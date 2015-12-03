from l2func_bind import *

class Add:
    def __init__(self, a, b):
        self.left = a
        self.right = b

class ScalarMult:
    def __init__(self, scalar, other):
        self.scalar = scalar
        self.other = other

class OpMult:
    def __init__(self, op, func):
        self.op = op
        self.func = func

# represent scalar        
class Sc:
    def __init__(self, value):
        self.value = value
    def __mul__(self, other):
        return ScalarMult(self.value, other)

def add(self, other):
    return Add(self, other)

def op_mult(self, func):
    return OpMult(self, func)

STO.__add__ = add
GTO.__add__ = add

D1.__add__ = add
D2.__add__ = add
Rm.__add__ = add

D1.__mul__ = op_mult
D2.__mul__ = op_mult
Rm.__mul__ = op_mult

cip_dict = {}
cip_dict[(STO, STO)] = cip_ss
cip_dict[(GTO, STO)] = cip_gs
cip_dict[(STO, GTO)] = cip_sg
cip_dict[(GTO, GTO)] = cip_gg

def get_cip(a, b):
    return cip_dict[(type(a), type(b))]

cip_op_dict = {}
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

def get_cip_op(a, o, b):

    return cip_op_dict[(type(a), type(o), type(b))]

def cip_op(a, o, b):
    """a and b is not Add neither ScalarMult"""
    if(isinstance(o, Add)):
        return  cip_op(a, o.left, b) + cip_op(a, o.right, b)
    if(isinstance(o, ScalarMult)):
        return o.scalar * cip_op(a, o.other, b)
    return get_cip_op(a, o, b)(a, o, b)

def cip(a, b):

    if(isinstance(a, Add)):
        return cip(a.left, b) + cip(a.right, b)
    if(isinstance(a, ScalarMult)):
        return a.scalar * cip(a.other, b)
    if(isinstance(b, Add)):
        return cip(a, b.left) + cip(a, b.right)
    if(isinstance(b, ScalarMult)):
        return b.scalar * cip(a, b.other)

    if(isinstance(a, OpMult)):
        if(isinstance(b, OpMult)):
            raise Exception("cip(OpMult, OpMult) is not supported now")
        return cip_op(a.func, a.op, b)
    
    if(isinstance(b, OpMult)):
        if(isinstance(a, OpMult)):
            raise Exception("cip(OpMult, OpMult) is not supported now")
        return cip_op(a, b.op, b.func)

    return get_cip(a, b)(a, b)


"""
def cip_ab(a, b):
    if(isinstance(a, Add)):
        return cip_ab(a.left, b) + cip_ab(a.right, b)
    if(isinstance(a, Mult)):
        return a.left * cip_ab(a.right, b)
    if(isinstance(b, Add)):
        return cip_ab(a, b.left) + cip_ab(a, b.right)
    if(isinstance(b, Mult)):
        return b.left * cip_ab(a, b.left)
    return get_cip(a, b)(a, b)

def cip_a_op_b(a, o, b):
    if(isinstance(a, Add)):
        return cip_a_op_b(a.left, o, b) + cip_a_op_b(a.right, o, b)
    if(isinstance(a, Mult)):
        return a.left * cip_a_op_b(a.right, o, b)
    if(isinstance(b, Add)):
        return cip_a_op_b(a, o, b.left) + cip_a_op_b(a, o, b.right)
    if(isinstance(b, Mult)):
        return b.left * cip_a_op_b(a, o, b.left)
    if(isinstance(o, Add)):
        return cip_a_op_b(a, o.left, b) + cip_a_op_b(a, o.right, b)
    if(isinstance(o, Mult)):
        return o.left * cip_a_op_b(a, o.right, b)
    return get_cip_op(a, o, b)(a, o, b)

def cip(a, b, c=None):
    if(c == None):
        return cip_ab(a, b)
    else:
        return cip_a_op_b(a, b, c)
"""


