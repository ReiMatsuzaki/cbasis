from sympy import *

r = Symbol("r")
z = Symbol("z")

for n in [1,2,3,4]:
    func = (r**n)*exp(-z*r*r)
    nterm = 1/sqrt(integrate(func*func, (r,0,oo), conds='none'))
    print n, "GTO:", nterm, func

""" output
1 GTO: 2*2**(3/4)/(pi**(1/4)*sqrt(z**(-3/2))) r*exp(-r**2*z)
2 GTO: 4*2**(3/4)*sqrt(3)/(3*pi**(1/4)*sqrt(z**(-5/2))) r**2*exp(-r**2*z)
3 GTO: 8*sqrt(15)*2**(3/4)/(15*pi**(1/4)*sqrt(z**(-7/2))) r**3*exp(-r**2*z)
4 GTO: 16*sqrt(105)*2**(3/4)/(105*pi**(1/4)*sqrt(z**(-9/2))) r**4*exp(-r**2*z)
"""    


    
