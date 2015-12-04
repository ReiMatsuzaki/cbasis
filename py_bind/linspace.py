# ==== Basic Component for Linear Space ====
# ---- Linear Algebra ----
class FuncAdd:

    def __init__(self, a, b):
        self.left = a
        self.right = b

    def __repr__(self):
        return "FuncAdd({0}, {1})".format(self.left, self.right)

    def __str__(self):
        return "{0} + {1}".format(self.left, self.right)

class ScalarFuncMult:
    def __init__(self, scalar, func):
        self.scalar = scalar
        self.func = func

    def __repr__(self):
        return "ScalarFuncMult({0}, {1})".format(self.scalar, self.func)

    def __str__(self):
        if(isinstance(self.func, FuncAdd)):
            return "{0}*({1})".format(self.scalar, self.func)
        else:
            return "{0}*{1}".format(self.scalar, self.func)

class OpAdd:
    def __init__(self, a, b):
        self.left = a
        self.right = b

    def __repr__(self):
        return "OpAdd({0}, {1})".format(self.left, self.right)

    def __str__(self):
        return "{0} + {1}".format(self.left, self.right)

class ScalarOpMult:
    def __init__(self, scalar, op):
        self.scalar = scalar
        self.op = op

    def __repr__(self):
        return "ScalarOpMult({0}, {1})".format(self.scalar, self.op)

    def __str__(self):
        if(isinstance(self.op, OpAdd)):
            return "{0}*({1})".format(self.scalar, self.op)
        else:
            return "{0}*{1}".format(self.scalar, self.op)

class OpMult:
    def __init__(self, op, func):
        self.op = op
        self.func = func

    def __repr__(self):
        return "OpMult({0}, {1})".format(self.scalar, self.func)

    def __str__(self):
        return "{0}[{1}]".format(self.scalar, self.func)

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
    return self.left.at(x) + self.right.at(x)

FuncAdd.at = func_add_at
def func_mult_scalar_at(self, x):
    return self.scalar * self.func.at(x)

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

"""
def cip_op(a, o, b):
    if(isinstance(o, OpAdd)):
        return  cip_op(a, o.left, b) + cip_op(a, o.right, b)
    if(isinstance(o, ScalarOpMult)):
        return o.scalar * cip_op(a, o.op, b)
    return get_cip_op(a, o, b)(a, o, b)
"""

def cip_impl(a, o, b):

    if(isinstance(a, FuncAdd)):
        return cip_impl(a.left, o, b) + cip_impl(a.right, o, b)
    if(isinstance(a, ScalarFuncMult)):
        return a.scalar * cip_impl(a.func, o, b)
    if(isinstance(o, OpAdd)):
        return cip_impl(a, o.left, b) + cip_impl(a, o.right, b)
    if(isinstance(o, ScalarOpMult)):
        return o.scalar * cip_impl(a, o.op, b)
    if(isinstance(b, FuncAdd)):
        return cip_impl(a, o, b.left) + cip_impl(a, o, b.right)
    if(isinstance(b, ScalarFuncMult)):
        return b.scalar * cip_impl(a, o, b.func)

    if(isinstance(o, OpId)):
        return get_cip(a, b)(a, b)
    else:
        return get_cip_op(a, o, b)(a, o, b)
    

def cip(a, b, c = None):
    
    if(c == None):
        """
        if(isinstance(b, OpMult) and isinstance(a, OpMult)):
            raise Exception("a:OpMult and b:OpMult is not supported")
        if(isinstance(b, OpMult)):
            return cip_impl(a, b.op, b.func)
        if(isinstance(a, OpMult)):
            return cip_impl(a.func, a.op, b)
        """
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

