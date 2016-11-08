import sys
from r1basis import *
import datetime
import pandas as pd

## ==== Utils ====
def print_timestamp(label, out):
    t = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    out.write("[{0}] {1}\n".format(label.rjust(10, ' '), t))
    
def switch(k0, dic, use_exception=True, default=None):
    for (key, val) in dic.items():
        if(k0 == key):
            return val
    if(use_exception):
        raise(Exception("key not found"))
    return default

## ==== derivative basis ====
def dict_set_default(dic, key, val):
    if(key not in dic):
        dic[key] = val
    
def one_deriv(us, opt_list):
    """ compute derivative basis for orbital exponent.
    see support/normalized_exp.py

    us : STOs or GTOs
    .         input basis set
    """
    if(not us.has_coef_all()):
        raise Exception("coef is not set")
    
    if(us.size() != len(opt_list)):
        raise Exception("opt_list : invalid size")

    if(us.exp_power() != 1 and us.exp_power() != 2):
        raise Exception("invalid input. Input type = {0}".format(type(us)))

    dus = STOs() if us.exp_power() == 1 else GTOs()
    
    for i in range(us.size()):
        if(us.is_prim(i) and us.is_normal(i) and opt_list[i]):
            ui = us.basis(i)
            ci = ui.c(0)
            ni = ui.n(0)
            zi = ui.z(0)
            if(us.exp_power() == 1):
                dui = LC_STOs()
                dui.add(-ci,           ni+1,   zi)
                dui.add(ci * (ni+0.5)/zi, ni, zi)
            else:
                dui = LC_GTOs()
                dui.add(-ci, ni+2, zi)
                dui.add(ci * (2*ni+1)/(4.0*zi), ni,   zi)
            dus.add_not_normal(dui)            

    dus.setup()
    return dus

    """
    if(us.exp_power() == 1):
        dus = STOs()
        for i in range(us.size()):
            ui = us.basis(i)
            ci = ui.c(0)
            ni = ui.n(0)
            zi = ui.z(0)

            dui = LC_STOs()
            dui.add(-ci, ni+1,   zi)
            dci = ci * (ni+0.5)/zi;
            dui.add(dci, ni, zi)
            
            dus.add_not_normal(dui)

        dus.setup()
        return dus
            
    elif(us.exp_power() == 2):
        dus = GTOs()
        for i in range(us.size()):
            ui = us.basis(i)
            ci = ui.c(0)
            ni = ui.n(0)
            zi = ui.z(0)
            
            dui = LC_GTOs()
            dci = ci * (2*ni+1)/(4.0*zi)
            dui.add(dci, ni,   zi)
            dui.add(-ci, ni+2, zi)

            dus.add_not_normal(dui)

        dus.setup()
        return dus
"""

def two_deriv(us, opt_list):
    """ compute second derivative basis for orbital exponent.
    see support/normalized_exp.py

    us : STOs or GTOs
    .         input basis set
    """
    if(not us.has_coef_all()):
        raise Exception("coef is not set")
    
    if(us.size() != len(opt_list)):
        raise Exception("opt_list : invalid size")

    if(us.exp_power() != 1 and us.exp_power() != 2):
        raise Exception("invalid input. Input type = {0}".format(type(us)))

    dus = STOs() if us.exp_power() == 1 else GTOs()

    for i in range(us.size()):
        if(us.is_prim(i) and us.is_normal(i) and opt_list[i]):
            ui = us.basis(i)
            ci = ui.c(0)
            ni = ui.n(0)
            zi = ui.z(0)
            if(us.exp_power() == 1):
                dui = LC_STOs()
                dui.add(ci,                         ni+2, zi)                
                dui.add(ci*(-2.0)*(ni+0.5)/zi,      ni+1, zi)
                dui.add(ci*(ni*ni-1.0/4.0)/(zi*zi), ni,   zi)
            else:
                dui = LC_GTOs()
                dui.add(ci,                             ni+4, zi)
                dui.add(ci*(-2)*(2*ni+1)/(4*zi),        ni+2, zi)
                dui.add(ci*(4*ni*ni-4*ni-3)/(16*zi*zi), ni,   zi)
            dus.add_not_normal(dui)            

    dus.setup()
    return dus            

    """    
    if(us.exp_power() == 1):
        ddus = STOs()
        for i in range(us.size()):
            ui = us.basis(i)
            ci = ui.c(0)
            ni = ui.n(0)
            zi = ui.z(0)

            dui = LC_STOs()
            dui.add(ci, ni+2,  zi)
            dci = ci * (-2.0) * (ni+0.5)/zi
            dui.add(dci, ni+1, zi)
            ddci= ci * (ni*ni - 1.0/4.0) / (zi * zi)
            dui.add(ddci, ni, zi)

            dus.add_not_normal(dui)

        dus.setup()
        return dus
            
    elif(us.exp_power() == 2):
        dus = GTOs()
        for i in range(us.size()):
            ui = us.basis(i)
            ci = ui.c(0)
            ni = ui.n(0)
            zi = ui.z(0)

            dui = LC_GTOs()
            dui.add(ci, ni+4,  zi)
            dci = ci * (-2) * (2*ni+1)/(4*zi)
            dui.add(dci, ni+2, zi)
            ddci= ci * (4*ni*ni - 4*ni -3) / (16 * zi * zi)
            dui.add(ddci, ni, zi)
            
            dus.add_not_normal(dui)

        dus.setup()
        return dus
    """

## ==== hydrogen atom photoionization ====
## ==== H photoionization ====
class H_Photoionization():
    def __init__(self, channel, dipole):

        (n0, l0, l1) = switch(channel,
                              {"1s->kp": (1, 0, 1),
                               "2p->ks": (2, 1, 0),
                               "2p->kd": (2, 1, 2),
                               "3d->kp": (3, 2, 1),
                               "3d->kf": (3, 2, 3)})
        
        self.E0 = -1.0 / (2.0 * n0 * n0)
        self.L1 = l1

        if(dipole == "length"):
            self.driv = switch(channel,
                               {"1s->kp": LC_STOs().add(2.0, 2, 1.0)})
        elif(dipole == "velocity"):
            self.driv = switch(channel,
                               {"1s->kp": LC_STOs().add(2.0, 1, 1.0)})
        else:
            raise(Exception("invalid dipole"))


    def h_mat(self, a, b):
        d2 = calc_d2_mat(a, b)
        h = -0.5 * d2
        if(self.L1 != 0):
            lam = self.L1 * (self.L1 + 1)
            r2 = calc_rm_mat(a, -2, b)
            h += 0.5 * lam * r2
            
        r1 = calc_rm_mat(a, -1, b)
        h += -r1
        return h

    def l_mat(self, w, a, b):
        ene = w + self.E0
        return calc_rm_mat(a, 0, b) * ene - self.h_mat(a, b)

    def s_mat(self, a, b):
        return calc_rm_mat(a, 0, b)
    
    def dip_vec(self, a):
        return calc_vec(a, self.driv)
    
## ==== optimization utility ====
def get_opt_index(opt_list):
    """
    [False, False]       -> []
    [True, True]         -> [0, 1]
    [True, False, True]  -> [0, 2]
    [False, False, True] -> [2]
    """
    
    num = len(opt_list)
    num_d = len([1 for opt_q in opt_list if opt_q])
    opt_index = range(num_d)
    i_d = 0
    for i in range(num):
        if opt_list[i]:
            opt_index[i_d] = i
            i_d = i_d + 1
    return opt_index

def update_basis(base_us, opt_index, zs):

    """
    Inputs
    -------
    base_us : STOs or GTOs
    opt_index : [int]
    .        obtained by get_opt_index
    zs      : current value of orbital exponents for optimizing basis
    """

    num   = base_us.size()
    num_d = len(zs)

    us = base_us.clone()
    for i_d in range(num_d):
        i = opt_index[i_d]
        bi = base_us.basis(i)
        bi.set_z(0, zs[i_d])
        us.replace(i, bi)
        us.setup()

    return us

## ==== Variables transform ====
class VarTrans():
    """
    Given the relation the sequence 
    .       {x_i | i = 1,...,Nx}, 
    and functions 
    .       {y_k(x) | k = 1,...,Ny}
    compute gradient/Hessian etc.

    dF/dyk = sum_i dxi/dyk dF/dxi
    d2F/(dyk,dyl) = sum_i d2xi /(dyk,dyl) dF/dxi + sum_ij dxi/dyk dxj/dyl H_ij

    """
    
    def __init__(self, Nx, Ny):
        self.Nx = Nx
        self.Ny = Ny

    def yks(self, xis):
        pass

    def xis(self, yks):
        pass
        
    def dF_dyk(self, dF_dyk, xis):
        return self.xis_yks(xis).transpose() * gis

    def d2F_dykdyl(self, dF_dxi, d2F_dxidxj, xis):
        raise(Exception("not implemented"))

    def dxi_dyk(self, xis):
        raise(Exception("not implemented"))
        
class VarTransComb(VarTrans):
    def __init__(self, xidx_list, var_list):
        self.xidx_list = xidx_list
        self.var_list  = var_list
        self.Nx = sum([var.Nx for var in var_list])
        
        self.yidx_list = []
        self.Ny = 0
        n = 0
        for var in var_list:
            self.yidx_list.append(range(n, n+var.Ny))
            n = n + var.Ny
            self.Ny = self.Ny + var.Ny

    def xis(self, yks):
        xis_res = VectorXc.Zero(self.Nx)
        for (xidx, yidx, var) in zip(self.xidx_list,
                                     self.yidx_list,
                                     self.var_list):
            yks_local = VectorXc([yks[i] for i in yidx])
            xis_local = var.xis(yks_local)
            for (i, x) in zip(xidx, xis_local):
                xis_res[i] = x
        return xis_res
    def yks(self, xis):
        yks_res = VectorXc.Zero(self.Ny)
        for (xidx, yidx, var) in zip(self.xidx_list,
                                     self.yidx_list,
                                     self.var_list):
            _xis = VectorXc([xis[i] for i in xidx])
            _yks = var.yks(_xis)
            for (i, y) in zip(yidx, _yks):
                yks_res[i] = y
        return yks_res
        
    def dF_dyk(self, dF_xi, xis):
        res = VectorXc.Zero(self.Ny)
        for (xidx, yidx, var) in zip(self.xidx_list,
                                     self.yidx_list,
                                     self.var_list):
            dF_xi_local = [dF_xi[i] for i in xidx]
            xis_local = [xis[i] for i in xidx]
            dF_dyk_local = var.dF_dyk(dF_xi_local, xis_local)
            for (i, dF) in zip(yidx, dF_dyk_local):
                res[i] = dF
        return res

    def d2F_dykdyl(self, dF_dxi, dF2_dxidxj, xis):
        res = MatrixXc.Zero(self.Ny, self.Ny)
        xyidx_var = zip(self.xidx_list, self.yidx_list, self.var_list)
        """
        for (xidx, yidx, var) in xyidx_var:
            _dF_dxi     = [dF_dxi[i] for i in xidx]
            _d2F_dxidxj = np.array([[dF2_dxidxj[i,j] for i in xidx] for j in xidx])
            _xis      = [xis[i]   for i in xidx]
            _dF_dyk= var.dF_dyk(_dF_dxi, _xis)
            _d2F_dykdyl= var.d2F_dykdyl(_dF_dxi, _d2F_dxidxj, _xis)
        for (k, _k)  in zip(yidx, range(var.Ny)):
                for (l, _l)  in zip(yidx, range(var.Ny)):
                    res[k, l] = _d2F_dykdyl[_k, _l]
        """
        num = len(xyidx_var)
        for ii in range(num):
            (xidx, yidx, vari) = xyidx_var[ii]
            _dF_dxi     = [dF_dxi[i] for i in xidx]
            _d2F_dxidxj = np.array([[dF2_dxidxj[i,j] for i in xidx]
                                    for j in xidx])
            _xis      = [xis[i]   for i in xidx]
            _dF_dyk= vari.dF_dyk(_dF_dxi, _xis)
            _d2F_dykdyl= vari.d2F_dykdyl(_dF_dxi, _d2F_dxidxj, _xis)
            for (k, _k)  in zip(yidx, range(vari.Ny)):
                for (l, _l)  in zip(yidx, range(vari.Ny)):
                    res[k, l] = _d2F_dykdyl[_k, _l]
            
            for jj in range(num):
                if ii!=jj:
                    (xjdx, yjdx, varj) = xyidx_var[jj]
                    _xjs = [xis[j] for j in xjdx]
                    _dxi_dyk = vari.dxi_dyk(_xis)
                    _dxj_dyl = varj.dxi_dyk(_xjs)
                    _d2F_dxidxj = np.array([[dF2_dxidxj[i,j] for j in xjdx]
                                            for i in xidx])
                    _d2F_dykdyl = _dxi_dyk.transpose() * _d2F_dxidxj * _dxj_dyl
                    """
                    print ii, jj
                    print "xidx,xjdx", xidx, xjdx
                    print "xis,xjs",   _xis, _xjs
                    print "dxi/dyk: ", np.array(_dxi_dyk)
                    print "dxj/dyl: ", np.array(_dxj_dyl)
                    print "dF2/ij : ", np.array(_d2F_dxidxj)
                    """
                    for (k, _k)  in zip(yidx, range(vari.Ny)):
                        for (l, _l)  in zip(yjdx, range(varj.Ny)):
                            res[k, l] = _d2F_dykdyl[_k, _l]            
                
        return res

class VarTransShift(VarTrans):
    """
    xi = a0k + yk[0]
    dF/dy0 = sum_i dxi/dy0 dF/dxi = sum_i dF/dxi
    dF/(dy0,dy0) =  sum_ij H_ij
    """
    def __init__(self, a0s):
        self.a0s = a0s
        self.Nx = len(a0s)
        self.Ny = 1

    def xis(self, yks):
        return VectorXc([a0+yks[0] for a0 in self.a0s])

    def yks(self, xis):
        return VectorXc([-self.a0s[0]+xis[0]])
    def xis(self, yks):
        return VectorXc([a0 + yks[0] for a0 in self.a0s])
    def dF_dyk(self, dF_dxi, xis):
        return VectorXc([sum(dF_dxi)])

    def d2F_dykdyl(self, dF_dxi, dF2_dxidxj, xis):
        ele = sum(sum(dF2_dxidxj))
        return MatrixXc([[ele]])

    def dxi_dyk(self, xis):
        res = MatrixXc.Zero(self.Nx, self.Ny)
        for i in range(self.Nx):
            res[i, 0] = 1.0
        return res
    
class VarTransMono(VarTrans):
    """
    yi(x) = f(xi)
    dF/dy = dxi/dyi dF/dxi
    dF2/dy2 = d2x/dy2 dF/dxi + (dxi/dyi)(dxi/dyi) dF2/dxi2
    """
    def __init__(self):
        self.Nx = 1
        self.Nx = 1

    def yks(self, xis):
        return VectorXc([self.y(xis[0])])

    def xis(self, yks):
        return VectorXc([self.x(yks[0])])    
        
    def dF_dyk(self, dF_dxi, xis):
        return VectorXc([self.dF_dy(dF_dxi[0], xis[0])])
    
    def dF_dy(self, dF_dx, x):
        return self.dx_dy(x) * dF_dx    

    def d2F_dykdyl(self, dF_dxi, d2F_dxidxj, xis):
        return MatrixXc([[self.d2F_dy2(dF_dxi[0], d2F_dxidxj[0,0], xis[0])]])
        
    def d2F_dy2(self, dF_dx, dF2_dx2, x):
        return self.dx2_dy2(x) * dF_dx + self.dx_dy(x)**2 * dF2_dx2
    
    def y(self, x):
        pass

    def x(self, y):
        pass

    def dx_dy(self, x):
        pass

    def dx2_dy2(self, x):
        pass


    def dxi_dyk(self, xis):
        return MatrixXc([[self.dx_dy(xis[0])]])
    
class VarTransId(VarTransMono):
    def y(self, x):
        return x

    def x(self, y):
        return y

    def dx_dy(self):
        return 1.0
    
class VarTransLog(VarTransMono):
    """
    yk = log(xi); xi = exp(yk)
    dxi/dyk = exp(yk)
    """
    def __init__(self):
        self.Nx = 1
        self.Ny = 1
    
    def y(self, x):
        return np.log(x)

    def x(self, y):
        return np.exp(y)
        
    def dx_dy(self, x):
        return x

    def dx2_dy2(self, x):
        return x

## ==== (Value,Grad,Hess) ====
def vgh_green(us, dus, ddus, opt_index, lmat, svec, rvec):
    num   = us.size()
    num_d = dus.size()

    L00 = lmat(us,  us) 
    L10 = lmat(dus,  us) 
    L11 = lmat(dus,  dus) 
    L20 = lmat(ddus,  us) 
    S0  = svec(us)
    S1  = svec(dus)
    S2  = svec(ddus)    
    R0  = rvec(us)
    R1  = rvec(dus)
    R2  = rvec(ddus)
            
    G = L00.inverse()

    val = tdot(S0, G*R0)
    grad = VectorXc.Zero(num_d)
    hess = MatrixXc.Zero(num_d, num_d)
    for i_d in range(num_d):
        i = opt_index[i_d]
        Si = VectorXc.Zero(num); Si[i] = S1[i_d]
        Ri = VectorXc.Zero(num); Ri[i] = R1[i_d]
        Li = MatrixXc.Zero(num,num);
        for ii in range(num):
            Li[i,ii] += L10[i_d, ii]
            Li[ii,i] += L10[i_d, ii]
        g_i = tdot(Si, G*R0) + tdot(S0, G*Ri) - tdot(S0, G*Li*G*R0)
        grad[i_d] = g_i
            
        for j_d in range(num_d):
            j = opt_index[j_d]
            Sj = VectorXc.Zero(num); Sj[j] = S1[j_d]
            Rj = VectorXc.Zero(num); Rj[j] = R1[j_d]
            Lj = MatrixXc.Zero(num, num);
            for jj in range(num):
                Lj[j,jj] += L10[j_d, jj]
                Lj[jj,j] += L10[j_d, jj]
            Lij= MatrixXc.Zero(num, num)
            if i==j:
                for jj in range(num):
                    Lij[jj,j] += L20[j_d, jj]
                    Lij[j,jj] += L20[j_d, jj]
            Lij[i,j] += L11[i_d, j_d]
            Lij[i,j] += L11[j_d, i_d]
            h_ij = (tdot(Sj, G*Ri) - tdot(Sj, G*Li*G*R0) + tdot(Si,G*Rj)
                    -tdot(S0, G*Lj*G*Ri) + tdot(S0, G*Lj*G*Li*G*R0)
                    -tdot(Si, G*Lj*G*R0) - tdot(S0, G*Lij*G*R0)
                    +tdot(S0, G*Li*G*Lj*G*R0) - tdot(S0, G*Li*G*Rj))
            if i==j:
                Sij = VectorXc.Zero(num); Sij[i]=S2[i_d]
                Rij = VectorXc.Zero(num); Rij[i]=R2[i_d]
                h_ij += tdot(Sij, G*R0) + tdot(S0, G*Rij)
            hess[i_d,j_d] = h_ij
    return (val, grad, hess)    
    
def vgh_green_h_pi(h_pi, base_us, opt_list):
    """ Returns value, gradient and Hessian for discretized matrix element of 
    Green's function. 
    .       <R|G|S>

    Inputs
    ------
    h_pi: H_Photoionization
    base_us : STOs or GTOs
    R  : LC_STO or LC_GTO
    S  : LC_STO or LC_GTO

    Returns
    -------
    val   : scalar
    grad  : VectorXc
    hess  : MatrixXc
    """
    if(base_us.size() != len(opt_list)):
        raise(Exception("size mismatch"))

    num   = base_us.size()
    num_d = len([1 for opt_q in opt_list if opt_q])
    opt_index = get_opt_index(opt_list)
    
    def __func_of_w__(w):
        def __func_of_zs__(zs):
            us = update_basis(base_us, opt_index, zs)
            dus = one_deriv(us, opt_list)
            ddus= two_deriv(us, opt_list)
            l_mat = lambda a,b: h_pi.l_mat(w, a, b)
            return vgh_green(us, dus, ddus, opt_index, 
                             l_mat, h_pi.dip_vec, h_pi.dip_vec)
        return __func_of_zs__
    return __func_of_w__

## ==== optimization ====

def vgh_log(vgh_w_zs):
    """
    y_i = log(z_i)  =>  z_i = exp(yi)
    df/dyi        = dzi/dyi gi = exp(yi) g_i
    d2f/(dyi dyj) = d/dyj exp(yi) g_i = exp(yi)exp(yj) g_ij
    d2f/dyi2 = d/dyi exp(yi) g_i = exp(yi) g_i + exp(2yi) g_ii
    """
    
    def __func_of_w__(w):
        def __func_of_ys__(ys):
            zs = [np.exp(y) for y in ys]
            (val, grad, hess) = vgh_w_zs(w)(zs)
            n = len(ys)
            gnew = VectorXc.Zero(n)
            hnew = MatrixXc.Zero(n, n)
            for i in range(n):
                gnew[i]   = zs[i] * grad[i]
                hnew[i,i] = zs[i] * grad[i] + zs[i]*zs[i]*hess[i,i]
                for j in range(n):
                    if i!=j:
                        hnew[i,j] = zs[i]*zs[j]*hess[i,j]
            return (val, gnew, hnew)
        return __func_of_ys__
    return __func_of_w__
            
class OptRes:
    def __init__(self):
        self.x = []          # solution array
        self.success = False # optimization successed or not
        self.message = ""    # calculation message
        self.val     = 0.0   # function values
        self.grad    = []    # function gradient
        self.hess    = [[]]  # its Hessian
        self.nit     = 0     # number of iterations

    def __str__(self):
        a = """
x = {0}
success = {1}
message = {2}
val     = {3}
grad    = {4}
hess    = {5}
nit     = {6}
""".format(self.x, self.success, self.message, self.val, self.grad, self.hess, self.nit)
        return a
        

## ==== optimization routine ====
def newton(vgh, x0, tol=0.00001, maxit=100):
    """
    Compute stationary point by Newton method.
    
    Inputs
    ------
    vgh : [scalar] -> (scalar, [scalar], [[scalar]])
    .     lambda for calculating function values, gradient and hessians
    x0  : [scalar]
    .     initial guess
    tol : real
    .     tolerrance for termination.    
    maxit : int
    .     maximum number of iterations

    Returns
    -------
    res : OptRes
    """

    res = OptRes()
    
    res.x = VectorXc(x0)
    for res.nit in range(maxit):
        (res.val, res.grad, res.hess) = vgh(res.x)
        ave_g = sum([abs(g0) for g0 in res.grad])/len(res.grad)
        if(ave_g < tol):
            res.success = True
            break
        dx =  -res.hess.inverse() * res.grad
        res.x = res.x + dx
        
    return res

## ==== Interface ====
def opt_main_init(args):

    ## ---- set default -----
    args['basis_type']; args['basis_info']
    args['w0'];  args['target']
    dict_set_default(args, 'tol',   pow(10.0, -5))
    dict_set_default(args, 'maxit', 10)
    dict_set_default(args, 'outfile', 'stdout')
    dict_set_default(args, 'wf_outfile', None)
    dict_set_default(args, 'ws', [args['w0']])
    dict_set_default(args, 'print_level', 0)
    dict_set_default(args, 'use_log', False)

    ## ---- check input ----
    if('ws_outfile' in args):
        if('ws_rs' in args):
            raise(Exception("ws_rs is necessary when ws_outfile is set"))

    ## ---- print input ----
    if(args['outfile'] == 'stdout'):
        args['out'] = sys.stdout
    else:
        args['out'] = open(args['outfile'], 'w')
    out = args['out']
    out.write('\n')
    out.write(">>>>opt_green>>>>>\n")
    out.write("optimize the orbital expoent for matrix element of Greens's operator: \n")
    out.write("       alpha(w) = <S, (E0+w-H)R> = SL^{-1}R\n")
    print_timestamp('Init', out)
    out.write("basis_type: {0}\n".format(args['basis_type']))
    out.write("basis_info:\n")
    for basis in args['basis_info']:
        out.write("{0}\n".format(basis))
    out.write("w0: {0}\n".format(args['w0']))
    out.write("ws: {0}\n".format(args['ws']))
    out.write("tol: {0}\n".format(args['tol']))
    out.write("target: {0}\n".format(args['target']))
    print >> out, "out = {0}".format("stdout" if args['outfile']
                                     else args['outfile'])
    print >> out, "wf_outfile = {0}".format(args['wf_outfile'] if args['wf_outfile']
                                            else "No")    
    
def opt_main_h_pi(args):

    ## ---- Check ----
    out = args["out"]
    print_timestamp('H_PI', out)
    out.write("channel: {0}\n".format(args['channel']))
    out.write("dipole:  {0}\n".format(args['dipole']))
    
    h_pi = H_Photoionization(args['channel'], args['dipole'])
    opt_list = []
    base_us = switch(args["basis_type"], {'STO': STOs(), 'GTO': GTOs()})
    for (pn, zeta, cmd) in args["basis_info"]:
        opt_list.append(switch(cmd, {'o': True, 'f': False}))
        base_us.add(pn, zeta)
    base_us.setup()
    vgh_w_zs = vgh_green_h_pi(h_pi, base_us, opt_list)
    z0s = [base_us.basis(i).z(0)
           for (opt, i)
           in zip(opt_list, range(base_us.size())) if opt]

    opt_index = get_opt_index(opt_list)
    args['h_pi'] = h_pi
    args['opt_index'] = opt_index
    args['base_us'] = base_us
    args['vgh_w_zs'] = vgh_w_zs
    args['z0s'] = z0s

def opt_main_calc(args):

    out = args['out']
    print_timestamp('Calc', out)
    tol = args['tol']
    maxit = args['maxit']    
    if(args['use_log']):
        vgh_w_zs = vgh_log(args['vgh_w_zs'])
        zs = args['z0s']
        n = zs.size()
        z0s  = VectorXc.Zero(n)
        for i in range(n):
            z0s[i] = np.log(zs[i])
    else:
        vgh_w_zs = args['vgh_w_zs'] 
        z0s  = args['z0s']
    
    ## ---- w range ----
    eps_w = pow(10.0, -5)
    w0 = args['w0']
    ws = args['ws']
    ws_without_w0 = [w for w in ws if abs(w-w0) > eps_w]
    ws_minus = [w0] + sorted([w for w in ws_without_w0 if w < w0],
                             key=lambda x:-x)
    ws_plus = sorted([w for w in ws_without_w0 if w > w0],
                     key=lambda x:+x)

    ## ---- calculation ----
    w_res_list = []
    zs = z0s
    for w in ws_minus:
        res = newton(vgh_w_zs(w), zs, tol=tol, maxit=maxit)
        w_res_list.append((w, res))
        if(args['use_log']):
            zs = [np.exp(x) for x in res.x]
        else:
            zs = res.x
            
        if(args['print_level'] > 0):
            out.write("w: {0}\n".format(w))
            for (i, z) in zip(range(len(zs)), zs):
                out.write("x{0}: {1}\n".format(i, z))
        if(args['print_level'] > 1):
            for (i, g) in zip(range(len(res.grad)), res.grad):
                out.write("grad{0}: {1}\n".format(i, g))

    zs = w_res_list[0][1].x
    for w in ws_plus:
        res = newton(vgh_w_zs(w), zs, tol=tol, maxit=maxit)
        w_res_list.append((w, res))
        if(args['use_log']):
            zs = [np.exp(x) for x in res.x]
        else:
            zs = res.x

        if(args['print_level'] > 0):
            out.write("opt for {0}\n".format(w))
            out.write("x: {0}\n".format(zs))            
        if(args['print_level'] > 1):
            out.write("grad: {0}\n".format(res.grad))

    w_res_list.sort(key = lambda wr: wr[0])
    args["w_res_list"] = w_res_list

def opt_main_wf_h_pi(args):

    out = args['out']
    print_timestamp("wf_h_pi", out)
    
    h_pi = args["h_pi"]
    base_us = args['base_us']
    w0 = args["w0"]
    (w0, res0) = [(w, res) for (w, res) in args["w_res_list"]
                  if abs(w-w0) < 0.00001][0]

    print >> out, "write wave function at w = {0}".format(w0)

    us = update_basis(base_us, args["opt_index"], res0.x)
    L00 = h_pi.l_mat(w0, us, us)
    R0  = h_pi.dip_vec(us)
    G = L00.inverse()
    cs = G * R0

    rs = args['wf_rs']
    ys = np.array(us.at_r(rs, cs))
    
    df = pd.DataFrame([rs, ys.real, ys.imag]).T
    df.columns = ["r", "re_y", "im_y"]
    df.to_csv(args['wf_outfile'], index=False)
    
def opt_main_print(args):
    out = args['out']
    
    print_timestamp('Result', out)
    for (w, res) in args['w_res_list']:
        out.write('\n')
        out.write("w = {0}\n".format(w))
        out.write("convergenec = {0}\n".format("Yes" if res.success else "No"))
        for i in range(len(res.x)):
            out.write("zeta{0} = {1}\n".format(i, res.x[i]))
        out.write("alpha = {0}\n".format(res.val))
        
def opt_main_finalize(args):

    print_timestamp('Finalize', args["out"])
    
    if(args["out"] != sys.stdout):
        args["out"].close()
        
def opt_main(**args):
    """ optimize 
    basis_type : str : STOs or GTOs
    basis_info : [(pn :int, zeta :complex, opt_cmd :str|other)]
    .             represent basis functions and its optimization command
    .   opt_cmd can be choosen as follow:
    .      'o'  :  optimize its orbital exponent
    .      'f'  :  fix its orbital exponent
    w0 : energy at start of optimization
    ws : list of w for optimization
    tol : tollerance in Newton method.
    target : str : calculation target. ex: "h_pi"
    out : file or None : str => output file. None => standard output
    wf_out : file or None : wave function output

    -- only if target = h_pi --
    channel : str : "1s->kp" or "2p->ks" etc
    dipole  : str : "length" or "velocity" etc
    """

    opt_main_init(args)
    if(args['target'] == "h_pi"):
        opt_main_h_pi(args)
    else:
        raise(Exception("only target='h_pi' is supported"))        
    opt_main_calc(args)
    opt_main_print(args)
    if(args['wf_outfile']):
        if(args['target'] == "h_pi"):
            opt_main_wf_h_pi(args)
    opt_main_finalize(args)
