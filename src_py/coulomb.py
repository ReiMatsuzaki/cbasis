import sympy
import cmath

def coulomb_phase(L, eta):
    """
    Gives the Coulomb phase shift.
    -------
    \sigma_L = arg Gamma(L+1+i\eta)
    \eta = z1 z2 / k
    -------
    """
    return cmath.phase(sympy.gamma(1+L+1.0j*eta))
    
