# ==== Basic Component for Linear Space ====
# ---- Base ----
class BaseAdd:
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def expand(self, f):
        return f(self.left) + f(self.right)

    def __str__(self):
        return "{0} + {1}".format(self.left, self.right)

class BaseScalarMult:
    def __init__(self, scalar, base):
        self.scalar = scalar
        self.base = base

    def expand(self, f):
        return self.scalar * f(self.base)

    def __str__(self):
        if(isinstance(self.base, BaseAdd)):
            return "{0}*({1})".format(self.scalar, self.base)
        else:
            return "{0}*{1}".format(self.scalar, self.base)


# ---- Func/Op ----        
class FuncAdd(BaseAdd): 
    def __repr__(self):
        return "FuncAdd({0}, {1})".format(self.left, self.right)

class ScalarFuncMult(BaseScalarMult):
    def __repr__(self):
        return "ScalarFuncMult({0}, {1})".format(self.scalar, self.base)

class OpAdd(BaseAdd):
    def __repr__(self):
        return "OpAdd({0}, {1})".format(self.left, self.right)

class ScalarOpMult(BaseScalarMult):
    def __repr__(self):
        return "ScalarOpMult({0}, {1})".format(self.scalar, self.base)

class OpMult:
    def __init__(self, op, func):
        self.op = op
        self.func = func

    def __repr__(self):
        return "OpMult({0}, {1})".format(self.scalar, self.func)

    def __str__(self):
        return "{0}[{1}]".format(self.scalar, self.func)


# ---- identity operator ----
class OpId:
    def __init__(self):
        pass

    def __repr__(self):
        return "OpId()"

    def __str__(self):
        return "id"


# ==== +,-,* ====
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

# ==== set operand ====
def set_as_func(FuncType):
    FuncType.__add__ = func_add
    FuncType.__neg__ = func_neg
    FuncType.__sub__ = func_sub
    FuncType.__rmul__ = scalar_func_mult
    FuncType.__mul__ = scalar_func_mult

map(set_as_func, [FuncAdd, ScalarFuncMult])

def set_as_op(OpType):
    OpType.__add__ = op_add
    OpType.__call__ = op_func_mult
    OpType.__neg__ = op_neg
    OpType.__sub__ = op_sub
    OpType.__rmul__ = scalar_op_mult
    OpType.__mul__ = scalar_op_mult

map(set_as_op, [OpAdd, ScalarOpMult])

# ==== set at ====
def func_add_at(self, x):
    return self.expand(lambda y: y.at(x))

FuncAdd.at = func_add_at
def func_mult_scalar_at(self, x):
    return self.expand(lambda y: y.at(x))

ScalarFuncMult.at = func_mult_scalar_at
    
# ==== Inner Product ====
cip_dict = {}
cip_op_dict = {}

def get_cip(a, b):
    try:
        the_cip = cip_dict[(type(a), type(b))]
    except KeyError:
        msg = "a: {0}, {1}\n".format(a.__repr__(), type(a))
        msg+= "b: {0}, {1}\n".format(b.__repr__(), type(b))
        raise KeyError(msg)
    return the_cip

def get_cip_op(a, o, b):
    try:
        the_cip = cip_op_dict[(type(a), type(o), type(b))]
    except KeyError:
        msg = "a: {0}, {1}\n".format(a.__repr__(), type(a))
        msg+= "o: {0}, {1}\n".format(o.__repr__(), type(o))
        msg+= "b: {0}, {1}\n".format(b.__repr__(), type(b))
        raise KeyError(msg)
    return the_cip


def cip_impl(a, o, b):

    if(isinstance(a, FuncAdd)):
        return a.expand(lambda x:cip_impl(x, o, b))
    if(isinstance(a, ScalarFuncMult)):
        return a.expand(lambda x:cip_impl(x, o, b))
    if(isinstance(o, OpAdd)):
        return o.expand(lambda x:cip_impl(a, x, b))
    if(isinstance(o, ScalarOpMult)):
        return o.expand(lambda x:cip_impl(a, x, b))
    if(isinstance(b, FuncAdd)):
        return b.expand(lambda x:cip_impl(a, o, x))
    if(isinstance(b, ScalarFuncMult)):
        return b.expand(lambda x:cip_impl(a, o, x))

    if(isinstance(o, OpId)):
        return get_cip(a, b)(a, b)
    else:
        return get_cip_op(a, o, b)(a, o, b)
    

def cip(a, b, c = None):
    
    if(c == None):
        return cip_impl(a, OpId(), b)
    else:
        return cip_impl(a, b, c)


# ==== induced from inner product ====
def cnorm2(f):
    return cip(f, f)

def cnorm(f):
    return sqrt(cnorm2(f))

def cnormalize(f):
    c = cnorm(f)
    return (1.0/c) * f

def linear_combination(cs, fs):
    """
    gives linear combination of function set fs with coefficients cs

    Inputs
    ------
    cs : [scalar]
    fs : [function]

    Returns
    -------
    res : function
    """

    def one(cumsum, cf):
        (c, f) = cf
        cumsum = cumsum + c*f

    return reduce(one, zip(cs, fs))

