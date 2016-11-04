from r1basis import *

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

class H_Photoionization():
    def __init__(self, n0, L0, L1, dipole):
        self.n0 = n0
        self.L0 = L0
        self.L1 = L1
        self.dipole = dipole
        self.E0 = -1.0 / (2.0 * n0 * n0)
        
        if(abs(L0-L1) != 1):
            raise(Exception("invalid channel"))

        if(dipole != "length" and dipole != "velocity"):
            raise(Exception("dipole <- {length, velocity}"))

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

    def l_mat(self, a, b, w):
        return calc_rm_mat(a, 0, b) * (w + self.E0) - self.h_mat(a, b)

    def dip_vec(self, a):
        if(self.dipole == "length"):
            if(self.n0 == 1 and self.L0 == 0):
                driv = LC_STOs().add(2.0, 2, 1.0)
            else:
                raise(Exception("not implemented yet"))
        else:
            if(self.n0 == 1 and self.L0 == 0):
                driv = LC_STOs().add(2.0, 1, 1.0)
            else:
                raise(Exception("not implemented yet"))
        return calc_vec(a, driv)


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

def vgh_green(h_pi, w, us, opt_list):
    """ Returns value, gradient and Hessian for discretized matrix element of 
    Green's function. 
    .       <R|G|S>

    Inputs
    ------
    h_pi: H_Photoionization
    us : STOs or GTOs
    R  : LC_STO or LC_GTO
    S  : LC_STO or LC_GTO

    Returns
    -------
    val   : scalar
    grad  : VectorXc
    hess  : MatrixXc
    """
    if(us.size() != len(opt_list)):
        raise(Exception("size mismatch"))

    num   = us.size()
    num_d = len([1 for opt_q in opt_list if opt_q])    
    opt_index = get_opt_index(opt_list)

    dus = one_deriv(us, opt_list)
    ddus= two_deriv(us, opt_list)
    L00 = h_pi.l_mat(us,  us, w)
    L10 = h_pi.l_mat(dus, us, w)
    L11 = h_pi.l_mat(dus, dus,w)
    L20 = h_pi.l_mat(ddus,us, w)
    R0  = h_pi.dip_vec(us)
    R1  = h_pi.dip_vec(dus)
    R2  = h_pi.dip_vec(ddus)
    S0  = h_pi.dip_vec(us)
    S1  = h_pi.dip_vec(dus)
    S2  = h_pi.dip_vec(ddus)

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
    
    
class OptRes:
    def __init__(self):
        self.x = []          # solution array
        self.success = False # optimization successed or not
        self.message = ""    # calculation message
        self.val     = 0.0   # function values
        self.grad    = []    # function gradient
        self.hess    = [[]]  # its Hessian
        self.nit     = 0     # number of iterations      
        
def newton(vgh, z0, tol, maxiter=100):
    """
    Compute stationary point by Newton method.
    
    Inputs
    ------
    vgh : [scalar] -> (scalar, [scalar], [[scalar]])
    .     lambda for calculating function values, gradient and hessians
    z0  : [scalar]
    .     initial guess
    tol : real
    .     tolerrance for termination.    
    maxit : int
    .     maximum number of iterations
    """

    print vgh(z0)
