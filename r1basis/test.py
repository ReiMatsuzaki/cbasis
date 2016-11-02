import sys
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eig, solve
from sympy import Symbol, oo, diff, integrate, exp, sqrt

from r1basis import *
from opt_driv import *
sys.path.append("../src_py/nnewton")
from nnewton import *

import unittest

class Test_driv_basis(unittest.TestCase):
    def setUp(self):
        pass

    def test_deriv(self):

        for (name, create) in [("STO", STOs), ("GTO", GTOs)]:
            n0 = 3
            c0 = [1.1]
            r0 = [2.5]
            z0 = [1.1]
            
            calc = lambda zs: create().add(3, zs[0]).setup().at_r(r0, c0)[0]

            stos = create()
            stos.add(n0, z0)
            stos.setup()
            
            dstos = one_deriv(stos)
            calc_dy = dstos.at_r(r0, c0)[0]
            ref_dy  = num_pd_c1(calc, z0, 0.0001, 0)        
            self.assertAlmostEqual(ref_dy, calc_dy)
            
            ddstos = two_deriv(stos)
            calc_ddy = ddstos.at_r(r0, c0)[0]
            ref_ddy  = num_pd2(calc, z0, 0.0001, 0, 0, method="c1")
            msg = "Basis={0}, second derivative\nref = {1}\ncalc= {2}\n".format(name, ref_ddy, calc_ddy)
            
            self.assertAlmostEqual(ref_ddy, calc_ddy, msg=msg)

            
class Test_first(unittest.TestCase):
    def setUp(self):
        pass

    def test_add(self):
        self.assertEqual(3, add(1,2))

class Test_r1_linear_comb(unittest.TestCase):
    def setUp(self):
        pass

    def test_size(self):

        for (m, fs) in zip([1, 2], [LC_STOs(), LC_GTOs()]): 
            fs.add(1.1, 2, 1.3-0.2j)
            fs.add(0.1, 1, 1.2-0.2j)
            fs.add(0.2, 3, 1.4-0.1j)
            self.assertEqual(3, fs.size())

    def test_getter(self):
        for (m, fs) in zip([1, 2], [LC_STOs(), LC_GTOs()]): 
            fs.add(1.1, 2, 1.3-0.2j)
            fs.add(0.1, 1, 1.2-0.2j)
            fs.add(0.2, 3, 1.4-0.1j)
            self.assertAlmostEqual(1.1, fs.c(0))
            self.assertAlmostEqual(1.2-0.2j, fs.z(1))
            self.assertEqual(3, fs.n(2))
        
    def test_conj(self):
        
        r = 1.4
        c0 = 1.1; n0 = 1; z0 = 0.35; 
        c1 = 1.2; n1 = 2; z1 = 0.3; 

        for (m, fs) in zip([1,2], [LC_STOs(), LC_GTOs()]): 
            fs.add(c0, n0, z0)
            fs.add(c1, n1, z1)
            y = fs.at_r([r])[0]
            self.assertAlmostEqual(y,
                                   c0 * r**n0 * np.exp(-z0**m)+
                                   c1 * r**n1 * np.exp(-z1**m))

    def test_clone(self):
        
        r = 1.4
        c0 = 1.1; n0 = 1; z0 = 0.35; 
        c1 = 1.2; n1 = 2; z1 = 0.3; 

        for (m, fs) in zip([1,2], [LC_STOs(), LC_GTOs()]): 
            fs.add(c0, n0, z0)
            fs.add(c1, n1, z1)
            
            fs2 = fs.clone()
            fs.add(1.1, 2, 3.3)
            self.assertEqual(2, fs2.size())

    def test_conj(self):
        
        r = 1.4
        c0 = 1.1; n0 = 1; z0 = 0.35; 
        c1 = 1.2; n1 = 2; z1 = 0.3; 

        for (m, fs) in zip([1,2], [LC_STOs(), LC_GTOs()]): 
            fs.add(c0, n0, z0)
            fs.add(c1, n1, z1)
            c_fs = fs.conj()
            y = fs.at_r([r])[0]
            cy= c_fs.at_r([r])[0]
            self.assertAlmostEqual(y.conjugate(), cy)

class Test_gto(unittest.TestCase):

    def setUp(self):
        gtos = GTOs()
        gtos.add(1, 1.1)
        gtos.add(2, 1.4)
        gtos.add(3, [1.2, 1.3-0.1j])

        g1 = LC_GTOs()
        g1.add(1.1, 5, 1.3)
        g1.add(0.1, 6, 1.4)
        gtos.add(g1)

        gtos.setup()
        self.gtos = gtos

        r = Symbol('r')
        self.r = r
        fs = [r**1*exp(-1.1*r**2),
              r**2*exp(-1.4*r**2),
              r**3*exp(-1.2*r**2),
              r**3*exp(-(1.3-0.1j)*r**2),
              1.1*r**5*exp(-1.3*r**2)+0.1*r**6*exp(-1.4*r**2)]
        self.fs = fs
        """
        fs = [lambda r: r     * np.exp(-1.1*r*r),
              lambda r: r*r   * np.exp(-1.4*r*r),
              lambda r: r*r*r * np.exp(-1.2*r*r),
              lambda r: r*r*r * np.exp(-(1.3-0.1j)*r*r),
              lambda r: 1.1*r**5*np.exp(-1.3*r*r) + 
              0.1*r**6*np.exp(-1.4*r*r)]
        self.fs = fs
        """
    
    def test_size(self):

        gtos = self.gtos
        self.assertEqual(5, gtos.size())
        self.assertFalse(gtos.only_prim())

        gtos = GTOs()
        gtos.add(2, 1.1)
        gtos.add(2, 1.1)
        gtos.add(3, [1.2, 1.3-0.1j])
        gtos.setup()
        self.assertEqual(4, gtos.size())
        self.assertTrue(gtos.only_prim())

    def test_accessor(self):

        gtos = self.gtos
        self.assertEqual(2, gtos.basis(1).n(0))
        self.assertEqual(6, gtos.basis(4).n(1))

    def test_clone(self):
        g = self.gtos.clone()
        self.gtos.add(2, 3.0)
        self.gtos.add(2, 3.0)
        self.gtos.add(2, 3.0)
        self.assertEqual(5, g.size())

    def test_conj(self):
        cg = self.gtos.conj()
        self.assertAlmostEqual(1.3+0.1j, cg.basis(3).z(0))
        
    def _test_raise(self):
        g = GTOs()
        g.add(2, 1.1)
        self.assertRaises(RuntimeError, calc_rm_mat(g, 0, g))

    def test_int_gto_lc(self):
        a = LC_GTOs()
        a.add(1.1, 2, 1.2)
        b = LC_GTOs()
        b.add(1.2, 3, 1.4)
        
        calc = int_lc(a, 2, b)
        ref,err  = quad(lambda r: (a.at_r([r])[0]*r*r*b.at_r([r])[0]).real, 0, 4.0)
        self.assertAlmostEqual(ref, calc)

    def test_at_r(self):
        n1 = 2; z1 = 1.0
        n2 = 2; z2 = 1.5
        
        g = GTOs()
        g.add(n1, z1)
        g.add(n2, z2)
        g.setup()
        
        c1 = g.basis(0).c(0)
        c2 = g.basis(1).c(0)
        
        d1 = 1.1
        d2 = 1.2
        
        r = 2.5

        self.assertAlmostEqual(d1*c1*(r**n1)*np.exp(-z1*r*r)+
                               d2*c2*(r**n2)*np.exp(-z2*r*r),
                               g.at_r([r], [d1, d2])[0])
        
    def test_hydrogen_atom(self):
        g =  GTOs()
        g.add(1, [2.5**n for n in range(-5,5)])
        g.setup()
        s = calc_rm_mat(g,0,g)
        h = -0.5 * calc_d2_mat(g, g) - calc_rm_mat(g, -1, g)
        (val,vec) =  eig(h, s)
        index = 5
        ene = val[index]
        c   = vec[:,index]
        self.assertAlmostEqual(-0.5, val[5], places=3)

        ## normalization
        nterm = 1.0/np.sqrt(np.dot(c, np.dot(s, c)))
        c1 = -nterm * c
        r0 = 2.5
        y_calc = g.at_r([r0], c1)[0]
        y_ref  = 2.0 * r0 * np.exp(-r0)
        self.assertAlmostEqual(y_ref, y_calc, places=4)

    def test_hydrogen_atom_p(self):
        g =  GTOs()
        g.add(2, [2.5**n for n in range(-5,5)])
        g.setup()
        s = calc_rm_mat(g, 0, g)
        h = -0.5 * calc_d2_mat(g,g) + calc_rm_mat(g,-2,g) - calc_rm_mat(g,-1,g)
        (val,vec) =  eig(h, s)
        self.assertAlmostEqual(-0.125, val[6], places=3)

class Test_sto(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_int_sto(self):
        z = 2.3
        for n in [0, 1, 3]:
            f = lambda r: r**n * np.exp(-z*r)
            nume, err = quad(f, 0, 15.0)
            self.assertAlmostEqual(nume, sto_int(n, z))

    def test_calc_vec(self):

        stos = STOs()
        stos.add(1, 1.1)
        stos.add(2, 1.2)
        stos.setup()

        sto1 = LC_STOs()
        sto1.add(2.0, 2, 1.0)
        
        vec = calc_vec(stos, sto1)
        self.assertEqual(2, len(vec))

    def test_hydrogen(self):
        s = STOs()
        s.add(2, 0.5)
        s.setup()
        
        h = -0.5 * calc_d2_mat(s,s) + calc_rm_mat(s,-2,s) - calc_rm_mat(s,-1,s)
        self.assertAlmostEqual(-0.125, h[0,0])

class Test_matele(unittest.TestCase):
    def setUp(self):
        pass

    def test_mat(self):
        s = STOs().add_not_normal(2.0, 3, 1.2-0.3j).setup()
        g = GTOs().add_not_normal(1.3, 2, 1.1-0.1j).setup()

        

        r2_s_lc = LC_STOs().add(2.0, 5, 1.2-0.3j)
        r2_g_lc = LC_GTOs().add(1.3, 4, 1.1-0.1j)
        
        s2s_calc = calc_rm_mat(s, 2, s)[0, 0]
        g2g_calc = calc_rm_mat(g, 2, g)[0, 0]
        sDs_calc = calc_d2_mat(s, s)[0, 0]
        gDg_calc = calc_d2_mat(g, g)[0, 0]

        s2s2 = calc_vec(s, r2_s_lc)[0]
        g2g2 = calc_vec(g, r2_g_lc)[0]

        ## see support/int_exp.py
        s2s_ref = -27.5296456055511 + 37.4411871137165j
        g2g_ref = 0.16651663387627 + 0.0546850960763247j
        sDs_ref = -0.526917955926652 - 1.46206690245985j
        gDg_ref = -0.395454606004856 - 0.0541117842324456j

        self.assertAlmostEqual(s2s_ref, s2s_calc)
        self.assertAlmostEqual(g2g_ref, g2g_calc)
        self.assertAlmostEqual(sDs_ref, sDs_calc)
        self.assertAlmostEqual(gDg_ref, gDg_calc)

        self.assertAlmostEqual(g2g_ref, g2g2)
        self.assertAlmostEqual(s2s_ref, s2s2)

class Test_driv(unittest.TestCase):

    def setUp(self):
        pass

    def test_alpha(self):

        ## from calc/stoh/l_5/res.d
        ss = STOs()

        ss.add(2, 0.9965751177-0.0013743026j)
        ss.add(2, 1.0030366528-0.2836728004j)
        ss.add(2, 0.8462928140-0.6952686244j)
        ss.add(2, 0.4818046345-1.0023929406j)
        ss.add(2, 0.1412093744-1.0662761427j)
        ss.setup()

        driv = LC_STOs()
        driv.add(2.0, 2, 1.0)
        
        ene = 0.5

        s = calc_rm_mat(ss, 0,  ss)
        d2= calc_d2_mat(ss,     ss)
        r2= calc_rm_mat(ss, -2, ss)
        r1= calc_rm_mat(ss, -1, ss)
        lmat = (  s * ene
                + d2* 0.5
                + r2* (-1.0)
                + r1)
        mvec = calc_vec(ss, driv)
        cs = solve(lmat, mvec)
        alpha = np.dot(cs, mvec)
        w = ene + 0.5
        ref = 1.88562800720386-0.362705406693342j
        self.assertAlmostEqual(3.0*ref, alpha)

        
"""
class Test_r1gtos(unittest.TestCase):
    def setUp(self):
        self.gtos = GTOs()
self.gtos.add(2, 1.1)
        self.gtos.add(3, [1.2, 1.3-0.1j])
        self.gtos.add_lc([(1.0, 1, 1.5), (1.1, 4, 1.4)])
        self.gtos.setup()

    def test_size(self):
        self.assertEqual(5, self.gtos.size_prim())
        self.assertEqual(4, self.gtos.size_basis())

    def test_param(self):
        self.assertEqual(2, self.gtos.n_prim(0))
        self.assertEqual(3, self.gtos.n_prim(1))
        self.assertEqual(3, self.gtos.n_prim(2))
        self.assertEqual(1, self.gtos.n_prim(3))
        self.assertEqual(4, self.gtos.n_prim(4))

        self.assertEqual(1.1, self.gtos.z_prim(0))
        self.assertEqual(1.2, self.gtos.z_prim(1))
        self.assertEqual(1.3-0.1j, self.gtos.z_prim(2))
        self.assertEqual(1.5, self.gtos.z_prim(3))
        self.assertEqual(1.4, self.gtos.z_prim(4))

    def test_set_conj(self):
        cg = self.gtos.conj()
        self.assertAlmostEqual(1.1, cg.z(0))
        self.assertAlmostEqual(1.3+0.1j, cg.z(2))
        (s, t, v) = cg.calc_mat_stv(self.gtos, 1)
        stos = R1STOs(); 
        stos.add(1.1, 2, 1.2)
        v1 = cg.calc_mat_sto(self.gtos, stos)
        print v1

    def test_matrix(self):
        (s,t,v) = self.gtos.calc_mat_stv(1)
        self.assertAlmostEqual(1.0, s[0, 0])
        self.assertAlmostEqual(1.0, s[1, 1])

        s2 = MatrixXc.Zero(1, 1)
        t2 = MatrixXc.Zero(1, 1)
        v2 = MatrixXc.Zero(1, 1)
        self.gtos.calc_mat_stv(1, s2, t2, v2)
        self.assertAlmostEqual(t2[0, 1], t[0, 1])

    def test_vector(self):
        stos = R1STOs()
        stos.add(1.1, 2, 1.2)
        stos.add(1.1-0.2j, 3, 1.3)
        m = self.gtos.calc_vec_sto(stos)

    def test_print(self):
        s = R1STO(1.0, 2, 1.1)
        print s

    def test_at_r(self):
        gs = R1GTOs()
        gs.add(2, 1.2)
        gs.normalize()
        print gs.at_r([1.3], [1.1])[0]
        print gs.deriv_at_r([1.3], [1.1])[0]
        print gs.deriv_2_at_r([1.3], [1.1])[0]


class Test_Opt(unittest.TestCase):
    def setUp(self):
        pass

    def test_opt(self):
        sto = R1STOs()
        sto.add(2, 1.1)
        driv = DrivSTO(sto)
        
        op = OpCoulomb(1, 0.5)
        
        zs = [0.2-0.1j, 0.8-0.5j]
        gs = R1GTOs()
        gs.add(2, zs)
        gs.normalize()

        target = OptAlpha(driv, op, gs)

        optimizer = OptNewton(100, 10.0**(-5), target, 0)
        res = optimizer.optimize(zs)
        
        self.assertTrue(res.conv_q)

    def test_opt_partial(self):
        sto = R1STOs(); sto.add(2, 1.0)
        driv = DrivSTO(sto);
        op   = OpCoulomb(1, 0.5)

        zs = [0.2-0.1j, 0.8-0.5j]
        gs = R1GTOs(); gs.add(2, zs); gs.normalize()
        target = OptAlphaPartial(driv, op, gs, [1])
        
        optimizer = OptNewton(100, 10.0**(-5), target, 0)
        res = optimizer.optimize([zs[1]])
        
        self.assertTrue(res.conv_q)

    def test_opt_alpha_shift(self):
        gs = R1GTOs()
        zs = [0.463925,
              1.202518,
              3.379649,
              10.6072,
              38.65163,
              173.5822,
              1170.498,
              0.16934112166516593,
              0.08989389391311804,
              0.055610873913491725,
              0.03776599632952126,
              0.02731159914174668,
              0.020665855224060142,
              0.016180602421004654,
              0.013011569667967734]
        gs.add(2, zs)
        gs.normalize()
        
        sto = R1STOs()
        sto.add(2.0, 2, 1.0)
        driv = DrivSTO(sto)
        
        op = OpCoulomb(1, 0.5)

        opt_idx = [7, 8, 9, 10, 11, 12, 13, 14];
        target = OptAlphaShift(driv, op, gs, opt_idx)
        optimizer = OptNewton(100, 10**(-5), target, 0)
        res = optimizer.optimize(-0.02j)
        
        self.assertTrue(res.conv_q)
        shift_ref = -0.00293368-0.0204361j
        self.assertAlmostEqual(shift_ref, res.zs[0], places=4);
        alpha_ref = -5.6568937518988989+1.0882823480377297j
        self.assertAlmostEqual(alpha_ref, res.val)

"""         
if __name__ == '__main__':
    unittest.main()
