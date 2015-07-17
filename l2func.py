from l2func_bind import *

def sto(c, n, z):
    res = STOs()
    res.add_one(1.0, STO(c, n, z))
    return res

"""
class STO:
    def sym_ip_with(other):
        return other.sym_ip_with_sto(self)
    def sym_ip_with_gto(gto):
        return sym_ip_sg(self, gto)
    def sym_ip_with_sto(sto):
        return sym_ip_ss(self, sto)

class GTO:
    def sym_ip_with(other):
        return other.sym_ip_with_sto(self)
    def sym_ip_with_gto(gto):
        return sym_ip_gg(self, gto)
    def sym_ip_with_sto(sto):
        return sym_ip_gs(self, sto)    

def sym_ip(a, b):
    return a.sym_ip_with(b)
"""



