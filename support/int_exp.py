from sympy import *

r = Symbol('r')
s = 2.0 * r**3 * exp(-(1.2-0.3j)*r)
g = 1.3 * r**2 * exp(-(1.1-0.1j)*r*r)
print "s2s_ref=", integrate(s*r*r*s, (r,0,oo)).evalf()
print "g2g_ref=", integrate(g*r*r*g, (r,0,oo)).evalf()
print "sDs_ref=", integrate(s*diff(s,r,2), (r,0,oo)).evalf()
print "gDg_ref=", integrate(g*diff(g,r,2), (r,0,oo)).evalf()

""" output
s2s_ref = -27.5296456055511 + 37.4411871137165*I
g2g_ref = 0.16651663387627 + 0.0546850960763247*I
sDs_ref = -0.526917955926652 - 1.46206690245985*I
gDg_ref = -0.395454606004856 - 0.0541117842324456*I
"""
