from l2func_bind import *

# ==== Basic Component for Linear Space ====
# ---- Linear Algebra ----
class FuncAdd:

    def __init__(self, a, b):
        self.left = a
        self.right = b

    def __repr__(self):
        return "{0} + {1}".format(self.left, self.right)

    def __str__(self):
        return "{0} + {1}".format(self.left, self.right)

class ScalarFuncMult:
    def __init__(self, scalar, func):
        self.scalar = scalar
        self.func = func

    def __repr__(self):
        if(isinstance(self.func, FuncAdd)):
            return "{0}*({1})".format(self.scalar, self.func)
        else:
            return "{0}*{1}".format(self.scalar, self.func)

    def __str__(self):
        return self.__repr__()

class OpAdd:
    def __init__(self, a, b):
        self.left = a
        self.right = b

    def __repr__(self):
        return "{0} + {1}".format(self.left, self.right)

    def __str__(self):
        return "{0} + {1}".format(self.left, self.right)


class ScalarOpMult:
    def __init__(self, scalar, op):
        self.scalar = scalar
        self.op = op

    def __repr__(self):
        if(isinstance(self.op, OpAdd)):
            return "{0}*({1})".format(self.scalar, self.op)
        else:
            return "{0}*{1}".format(self.scalar, self.op)

    def __str__(self):
        return self.__repr__()


class OpMult:

    def __init__(self, op, func):
        self.op = op
        self.func = func

    def __repr__(self):
        return "{0}({1})".format(self.scalar, self.func)

    def __str__(self):
        return "{0}[{1}]".format(self.scalar, self.func)


# ---- +,-,* ----
def func_add(self, other):
    return FuncAdd(self, other)

def scalar_func_mult(self, scalar):
    return ScalarFuncMult(scalar, self)

def func_neg(self):
    return ScalarFuncMult(-1.0, self)

def func_sub(self, other):
    return FuncAdd(self, func_neg(other))

def op_add(self, other):
    return OpAdd(self, other)

def scalar_op_mult(self, scalar):
    return ScalarOpMult(scalar, self)

def op_func_mult(self, func):
    return OpMult(self, func)

def op_neg(self):
    return ScalarFuncMult(-1.0, self)

def op_sub(self, other):
    return OpAdd(self, op_neg(other))


# ==== Set Func ====
# ---- general ----
def set_as_func(FuncType):
    FuncType.__add__ = func_add
    FuncType.__neg__ = func_neg
    FuncType.__sub__ = func_sub
    FuncType.__rmul__ = scalar_func_mult
    FuncType.__mul__ = scalar_func_mult

map(set_as_func, [GTO, STO, FuncAdd, ScalarFuncMult])

# ---- at ----
def func_add_at(self, x):
    return self.left.at(x) + self.right.at(x)

FuncAdd.at = func_add_at
def func_mult_scalar_at(self, x):
    return self.scalar * self.func.at(x)

ScalarFuncMult.at = func_mult_scalar_at

# ---- repr/str ----
def STO_repr(self):
    return "STO({0},{1},{2})".format(self.c, self.n, self.z)

def GTO_repr(self):
    return "GTO({0},{1},{2})".format(self.c, self.n, self.z)

STO.__repr__ = STO_repr
GTO.__repr__ = GTO_repr

# ==== Set Operator ====
# ---- general ----
def set_as_op(OpType):
    OpType.__add__ = op_add
    OpType.__call__ = op_func_mult
    OpType.__neg__ = op_neg
    OpType.__sub__ = op_sub
    OpType.__rmul__ = scalar_op_mult
    OpType.__mul__ = scalar_op_mult

map(set_as_op, [D1, D2, Rm, OpAdd, ScalarOpMult])

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
    
    if(isinstance(o, OpAdd)):
        return  cip_op(a, o.left, b) + cip_op(a, o.right, b)
    if(isinstance(o, ScalarOpMult)):
        return o.scalar * cip_op(a, o.op, b)
    return get_cip_op(a, o, b)(a, o, b)

def cip(a, b):

    if(isinstance(a, FuncAdd)):
        return cip(a.left, b) + cip(a.right, b)
    if(isinstance(a, ScalarFuncMult)):
        return a.scalar * cip(a.func, b)
    if(isinstance(b, FuncAdd)):
        return cip(a, b.left) + cip(a, b.right)
    if(isinstance(b, ScalarFuncMult)):
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


