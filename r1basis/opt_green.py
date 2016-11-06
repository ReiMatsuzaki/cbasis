from r1basis import *

## ==== Utils ====
def switch(k0, dic, use_exception=True, default=None):
    for (key, val) in dic.items():
        if(k0 == key):
            return val
    if(use_exception):
        raise(Exception("key not found"))
    return default

## ==== derivative basis ====
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
    
def vgh_green(us, dus, ddus, opt_index, ene, hmat, smat, svec, rvec):
    num   = us.size()
    num_d = dus.size()

    L00 = smat(us,  us) * ene - hmat(us, us)
    L10 = smat(dus,  us) * ene - hmat(dus, us)
    L11 = smat(dus,  dus) * ene - hmat(dus, dus)
    L20 = smat(ddus,  us) * ene - hmat(ddus, us)
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
#            print 'energy=', h_pi.E0, w, h_pi.E0+w
            return vgh_green(us, dus, ddus, opt_index, h_pi.E0+w,
                             h_pi.h_mat, h_pi.s_mat, h_pi.dip_vec, h_pi.dip_vec)
        return __func_of_zs__
    return __func_of_w__

## ==== optimization results ====
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
def newton(vgh, x0, tol=0.00001, maxiter=100):
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
    for res.nit in range(maxiter):
        (res.val, res.grad, res.hess) = vgh(res.x)
        ave_g = sum([abs(g0) for g0 in res.grad])/len(res.grad)
        if(ave_g < tol):
            res.success = True
            break
        dx =  -res.hess.inverse() * res.grad
        res.x = res.x + dx
        
    return res


