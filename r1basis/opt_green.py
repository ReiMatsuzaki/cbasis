from r1basis import *

def one_deriv(us):
    """ compute derivative basis for orbital exponent.
    see support/normalized_exp.py

    us : STOs or GTOs
    .         input basis set
    """

    if(not us.is_normal()):
        raise Exception("us must be normalized.")

    if(not us.only_prim()):
        raise Exception("only support non contracted basis.")

    if(us.exp_power() != 1 and us.exp_power() != 2):
        raise Exception("invalid input. Input type = {0}".format(type(us)))

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

def two_deriv(us):
    """ compute second derivative basis for orbital exponent.
    see support/normalized_exp.py

    us : STOs or GTOs
    .         input basis set
    """

    if(not us.is_normal()):
        raise Exception("us must be normalized.")

    if(not us.only_prim()):
        raise Exception("only support non contracted basis.")

    if(us.exp_power() != 1 and us.exp_power() != 2):
        raise Exception("invalid input. Input type = {0}".format(type(us)))

    if(us.exp_power() == 1):
        dus = STOs()
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
def vgh_green(us, R, S):
    """ Returns value, gradient and Hessian for discretized matrix element of 
    Green's function. 
    .       <R|G|S>

    Inputs
    ------
    us : STOs or GTOs
    R  : LC_STO or LC_GTO
    S  : LC_STO or LC_GTO

    Returns
    -------
    val   : scalar
    grad  : VectorXc
    hess  : MatrixXc
    """

    
    
    
