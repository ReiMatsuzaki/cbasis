import numpy as np
from r1basis import *
from scipy.integrate import quad
from scipy.linalg import eig, solve

import unittest

class Test_first(unittest.TestCase):
    def setUp(self):
        pass

    def test_add(self):
        print add(1,2)

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

class _Test_gto(unittest.TestCase):

    def setUp(self):
        gtos = GTOs()
        gtos.add(1, 1.1)
        gtos.add(2, 1.4)
        gtos.add(3, [1.2, 1.3-0.1j])

        g1 = LC_GTOs()
        g1.add(1.1, 5, 1.3)
        g1.add(0.1, 6, 1.4)
        gtos.add_lc(g1)

        gtos.setup()
        
        self.gtos = gtos
    
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
        
    def test_raise(self):
        g = GTOs()
        g.add(2, 1.1)
        self.assertRaises(RuntimeError, g.calc_rm_mat, 0)

    def test_int_gto(self):
        n = 2
        z = 2.3
        f = lambda r: r**n*np.exp(-z*r*r)
        nume, err = quad(f, 0, 4.0)
        self.assertAlmostEqual(nume, gto_int(n, z))

    def test_int_gto_lc(self):
        a = LC_GTOs()
        a.add(1.1, 2, 1.2)
        b = LC_GTOs()
        b.add(1.2, 3, 1.4)
        
        calc = int_lc(a, 2, b)
        ref,err  = quad(lambda r: (a.at_r([r])[0]*r*r*b.at_r([r])[0]).real, 0, 4.0)
        self.assertAlmostEqual(ref, calc)
        
    def test_rm_matrix(self):

        gs = self.gtos
        s  = gs.calc_rm_mat(0)
        r2 = gs.calc_rm_mat(2)
        self.assertAlmostEqual(1.0, s[0, 0])

        fs = [lambda r: r     * np.exp(-1.1*r*r),
              lambda r: r*r   * np.exp(-1.4*r*r),
              lambda r: r*r*r * np.exp(-1.2*r*r),
              lambda r: r*r*r * np.exp(-(1.3-0.1j)*r*r),
              lambda r: 1.1*r**5*np.exp(-1.3*r*r) + 
              0.1*r**6*np.exp(-1.4*r*r)]
        
        s00, err = quad(lambda r: fs[0](r)*fs[0](r), 0, 10.0)
        s11, err = quad(lambda r: fs[1](r)*fs[1](r), 0, 10.0)
        s44, err = quad(lambda r: fs[4](r)*fs[4](r), 0, 10.0)
        r2_01, err = quad(lambda r: fs[0](r)*fs[1](r)*r*r, 0, 10.0)
        r2_14, err = quad(lambda r: fs[1](r)*fs[4](r)*r*r, 0, 10.0)
        self.assertAlmostEqual(r2[0, 1], r2_01/np.sqrt(s00*s11))
        self.assertAlmostEqual(r2[1, 4], r2_14/np.sqrt(s11*s44))

    def test_d2_mat(self):
        gs = self.gtos
        d2 = gs.calc_d2_mat()
        dif= d2-d2.transpose()
        self.assertAlmostEqual(0.0, dif[0, 1])
        
    def test_hydrogen_atom(self):
        g =  GTOs()
        g.add(1, [2.5**n for n in range(-5,5)])
        g.setup()
        s = g.calc_rm_mat(0)
        h = -0.5 * g.calc_d2_mat() - g.calc_rm_mat(-1)
        (val,vec) =  eig(h, s)
        self.assertAlmostEqual(-0.5, val[5], places=3)
        print g

    def test_hydrogen_atom_p(self):
        g =  GTOs()
        g.add(2, [2.5**n for n in range(-5,5)])
        g.setup()
        s = g.calc_rm_mat(0)
        h = -0.5 * g.calc_d2_mat() + g.calc_rm_mat(-2) - g.calc_rm_mat(-1)
        (val,vec) =  eig(h, s)
        print val
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
        
        vec = stos.calc_vec(sto1)
        self.assertEqual(2, len(vec))
        

    def test_hydrogen(self):
        s = STOs()
        s.add(2, 0.5)
        s.setup()
        
        h = -0.5 * s.calc_d2_mat() + s.calc_rm_mat(-2) - s.calc_rm_mat(-1)
        self.assertAlmostEqual(-0.125, h[0,0])

class _Test_driv(unittest.TestCase):

    def setUp(self):
        pass

    def test_alpha(self):

        ## from calc/stoh/l_5/res.d
        ss = STOs()

        ss.add(2, 0.9965751177-0.0013743026j)
        ss.add(2, 1.0030366528  -0.2836728004j)
        ss.add(2, 0.8462928140  -0.6952686244j)
        ss.add(2, 0.4818046345  -1.0023929406j)
        ss.add(2, 0.1412093744  -1.0662761427j)
        ss.setup()

        driv = LC_STOs()
        driv.add(2.0, 2, 1.0)

        ene = 0.5

        s = ss.calc_rm_mat(0)
        d2= ss.calc_d2_mat()
        r2= ss.calc_rm_mat(-2)
        r1= ss.calc_rm_mat(-1)
        #print "s:",  s
        #print "d2:", d2
        lmat = (  s * ene
                + d2* 0.5
                + r2* (-1.0)
                + r1)
        #print "lmat:", lmat
        mvec = ss.calc_vec(driv)
        #print "mvec", mvec
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
