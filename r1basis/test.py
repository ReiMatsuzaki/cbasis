import sys
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eig, solve
from sympy import Symbol, oo, diff, integrate, exp, sqrt

from r1basis import *
from opt_green import *
sys.path.append("../src_py/nnewton")
from nnewton import *

import unittest
class Test_tmp(unittest.TestCase):
    def setUp(self):
        pass
    
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
        self.assertFalse(gtos.is_prim_all())

        gtos = GTOs()
        gtos.add(2, 1.1)
        gtos.add(2, 1.1)
        gtos.add(3, [1.2, 1.3-0.1j])
        gtos.setup()
        self.assertEqual(4, gtos.size())
        self.assertTrue(gtos.is_prim_all())

    def test_accessor(self):

        gtos = self.gtos
        self.assertEqual(2, gtos.basis(1).n(0))
        self.assertEqual(6, gtos.basis(4).n(1))
        basis_i = gtos.basis(0).clone()
        basis_i.set_z(0, 1.3)
        gtos.replace(0, basis_i)
        self.assertAlmostEqual(1.3, gtos.basis(0).z(0))
        
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

    def test_accessor(self):
        s2 = LC_STOs().add(1.2, 2, 1.1).add(1.1, 3, 1.4)
        fs = STOs().add(1, 2.0).add_not_normal(1.1, 2, 3.0).add(s2).setup()
        self.assertTrue(fs.is_prim(0))
        self.assertTrue(fs.is_prim(1))
        self.assertFalse(fs.is_prim(2))
        
        self.assertTrue(fs.has_coef(0))
        self.assertTrue(fs.has_coef(1))
        self.assertTrue(fs.has_coef(2))

        self.assertTrue(fs.is_normal(0))
        self.assertFalse(fs.is_normal(1))
        self.assertTrue(fs.is_normal(2))        

        fs = STOs().add(2, 1.1)
        self.assertFalse(fs.has_coef(0))
        
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

    def test_sto_gto(self):
        for n in [0,1,2,3,4,5]:
            s = STOs().add_not_normal(2.0, n, 1.2).setup()
            g = GTOs().add_not_normal(1.1, 0, 1.1).setup()
            ele = calc_rm_mat(s, 0, g)[0, 0]
            sg = lambda x: x**n*2.0*1.1*np.exp(-1.2*x-1.1*x*x)
            (ref, err) = quad(sg, 0, 10.0)
            self.assertAlmostEqual(ref, ele)
    
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

    def test_mat2(self):
        s = STOs().add_not_normal(2.2, 3, 1.1).setup()
        g = GTOs().add_not_normal(1.3, 2, 1.2).setup()
        
        s2s_calc = calc_rm_mat(s, 2, s)[0, 0]
        g2g_calc = calc_rm_mat(g, 2, g)[0, 0]
        s2g_calc = calc_rm_mat(s, 2, g)[0, 0]
        sDs_calc = calc_d2_mat(s, s)[0, 0]
        gDg_calc = calc_d2_mat(g, g)[0, 0]
        sDg_calc = calc_d2_mat(s, g)[0, 0]

        ## see support/int_exp.py
        s2s_ref= 161.644807242673
        g2g_ref= 0.131127436620057
        s2g_ref= 0.663645309086432
        sDs_ref= -3.38091660405710
        gDg_ref= -0.352470549634713
        sDg_ref= 0.208872645967760

        self.assertAlmostEqual(s2s_ref, s2s_calc)
        self.assertAlmostEqual(g2g_ref, g2g_calc)
        self.assertAlmostEqual(s2g_ref, s2g_calc)
        self.assertAlmostEqual(sDs_ref, sDs_calc)
        self.assertAlmostEqual(gDg_ref, gDg_calc)
        self.assertAlmostEqual(sDg_ref, sDg_calc)
        
class Test_green(unittest.TestCase):

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

    def test_deriv_one(self):

        for (name, create) in [("STO", STOs), ("GTO", GTOs)]:
            n0 = 3
            c0 = [1.1]
            r0 = [2.5]
            z0 = [1.1]
            
            calc = lambda zs: create().add(3, zs[0]).setup().at_r(r0, c0)[0]

            stos = create()
            stos.add(n0, z0)
            stos.setup()
            
            dstos = one_deriv(stos, [True])
            calc_dy = dstos.at_r(r0, c0)[0]
            ref_dy  = num_pd_c1(calc, z0, 0.0001, 0)        
            self.assertAlmostEqual(ref_dy, calc_dy)

            ddstos = two_deriv(stos, [True])
            calc_ddy = ddstos.at_r(r0, c0)[0]
            ref_ddy  = num_pd2(calc, z0, 0.0001, 0, 0, method="c1")
            msg = "Basis={0}, second derivative\nref = {1}\ncalc= {2}\n".format(name, ref_ddy, calc_ddy)
            
            self.assertAlmostEqual(ref_ddy, calc_ddy, msg=msg)

    def test_deriv_three(self):
        z0 = 0.52
        r0 = 2.6
        c0 = 1.2
        #basis= lambda z: GTOs().add(3, [1.1, z, 3.3]).setup()
        basis= lambda z: GTOs().add(3, 1.1).add(3, z).add(4, 3.3).setup()

        ds = one_deriv(basis(z0), [False, True, False])
        self.assertEqual(1, ds.size())

        calc = lambda z: basis(z).at_r(r0, [0,c0,0])
        ref_dy  = num_pd(calc, [z0], 0.0001, 0, method='c1')
        calc_dy = ds.at_r(r0, [c0])
        self.assertAlmostEqual(ref_dy, calc_dy)
        
    def test_grad_one(self):
        z0 = [1.3-0.2j]
        h_pi = H_Photoionization('1s->kp', "length")
        basis  = STOs().add(2, z0).setup()
        w = 1.0
        calc_green = lambda z: vgh_green_h_pi(h_pi, basis, [True])(w)(z)[0]
        (val, calc_grad, calc_hess) = vgh_green_h_pi(h_pi, basis, [True])(w)(z0)
        ref_grad = num_pd(calc_green,  z0, 0.0001, 0,    method='c1')
        ref_hess = num_pd2(calc_green, z0, 0.0001, 0, 0, method='c1')
        self.assertAlmostEqual(calc_grad[0],   ref_grad)
        self.assertAlmostEqual(calc_hess[0,0], ref_hess, places=6)

    def test_opt_index(self):
        self.assertEqual([],     get_opt_index([False, False]))
        self.assertEqual([0, 1], get_opt_index([True, True]))
        self.assertEqual([1, 3], get_opt_index([False, True, False, True]))
        self.assertEqual([0, 2], get_opt_index([True, False, True]))
        
    def test_gh_two(self):

        z0 = [1.3-0.2j, 0.5-0.9j]
        opt = [True, True]
        h_pi = H_Photoionization('1s->kp', "length")
        basis  = STOs().add(2, z0).setup()
        w = 1.0
        calc_green = lambda z: vgh_green_h_pi(h_pi, basis, opt)(w)(z)[0]
        
        (val, calc_grad, calc_hess) = vgh_green_h_pi(h_pi, basis, opt)(w)(z0)
        ref_grad = ngrad(calc_green, z0, 0.0001, method='c1')
        ref_hess = nhess(calc_green, z0, 0.0001, method='c1')
        for i in range(2):
            self.assertAlmostEqual(ref_grad[i], calc_grad[i])
            for j in range(2):
                self.assertAlmostEqual(ref_hess[i,j],
                                       calc_hess[i,j], places=6)

    def test_gh_two_of_three_easy(self):
        zs0 = [1.3-0.1j, 2.3-0.5j]
        z2 = 0.4-0.4j
        opt = [True, True, False]
        h_pi = H_Photoionization('1s->kp', "velocity")
        basis = STOs().add(2, 1.0).add(2,1.0).add(2,z2).setup()
        w = 0.9
        
        calc_green = lambda z: vgh_green_h_pi(h_pi, basis, opt)(w)(z)[0]

        (v, calc_g, calc_h) = vgh_green_h_pi(h_pi, basis, opt)(w)(zs0)
        ref_g = ngrad(calc_green, zs0, 0.0001, method='c1')
        ref_h = nhess(calc_green, zs0, 0.0001, method='c1')
        for i in range(2):
            self.assertAlmostEqual(ref_g[i], calc_g[i], places=5)
            for j in range(2):
                msg = "ref= {0}\ncalc={1}\n(i,j)=({2},{3})".format(ref_h[i,j],calc_h[i,j],i,j)
                self.assertAlmostEqual(ref_h[i,j], calc_h[i,j], places=5, msg=msg)

    def test_gh_one_of_two(self):
        zs = [1.3-0.1j, 2.3-0.5j]
        h_pi = H_Photoionization('1s->kp', "velocity")
        basis = STOs().add(2, zs).setup()
        w = 1.1

        (v, calc_g, calc_h) = vgh_green_h_pi(h_pi, basis, [False, True])(w)([zs[1]])
        (v, full_g, full_h) = vgh_green_h_pi(h_pi, basis, [True,  True])(w)(zs)

        self.assertAlmostEqual(full_g[1], calc_g[0], places=5)
        self.assertAlmostEqual(full_h[1,1], calc_h[0,0], places=5)
                
    def test_grad_two_of_three(self):
        opt_zs0 = [1.3-0.1j, 2.3-0.5j]
        z1 = 0.4-0.4j
        zs_all = [opt_zs0[0], z1, opt_zs0[1]]
        opt = [True,      False,  True]
        h_pi = H_Photoionization('1s->kp', "velocity")
        basis =  STOs().add(2, 1.0).add(2,z1).add(2, 1.0).setup()
        w = 1.0
        calc_green = lambda z: vgh_green_h_pi(h_pi, basis, opt)(w)(z)[0]

        (v, calc_g, calc_h) = vgh_green_h_pi(h_pi, basis, opt)(w)(opt_zs0)
        (v, full_g, full_h) = vgh_green_h_pi(h_pi, basis, [True,True,True])(w)(zs_all)
        ref_g = ngrad(calc_green, opt_zs0, 0.0001, method='c1')
        ref_h = nhess(calc_green, opt_zs0, 0.0001, method='c1')

        self.assertAlmostEqual(full_h[0,0], calc_h[0,0], places=5)
        self.assertAlmostEqual(full_h[0,2], calc_h[0,1], places=5)
        self.assertAlmostEqual(full_h[2,0], calc_h[1,0], places=5)
        self.assertAlmostEqual(full_h[2,2], calc_h[1,1], places=5)
        
        for i in range(2):
            self.assertAlmostEqual(ref_g[i], calc_g[i], places=5)            
            for j in range(2):
                msg = "ref= {0}\ncalc={1}\n(i,j)=({2},{3})".format(ref_h[i,j],calc_h[i,j],i,j)
                self.assertAlmostEqual(ref_h[i,j], calc_h[i,j], places=5, msg=msg)

        
    def test_opt_one(self):
        h_pi = H_Photoionization('1s->kp', "velocity")
        z0s = [0.6-0.6j]
        opt = [True for z in z0s]
        basis = STOs().add(2, z0s).setup()
        w = 1.0
        
        res = newton(vgh_green_h_pi(h_pi, basis, opt)(w), z0s)

        ## see calc/stoh/v_1/res.d
        self.assertTrue(res.success)
        self.assertAlmostEqual(1.0255886472-0.6955918398j, res.x[0])
        self.assertAlmostEqual((0.361600808054165-0.371221793708147j)*3, res.val)

    def test_opt_three(self):
        h_pi = H_Photoionization('1s->kp', "length")
        zs_opt = [0.9797019427  -0.0125136315j,
                  0.8771210224  -0.6400667900j,
                  0.3008012645  -1.0095895471j]
        z0s = [z0 + 0.01 for z0 in zs_opt]
        opt = [True for z in z0s]
        basis = STOs().add(2, z0s).setup()
        w = 0.9
        res = newton(vgh_green_h_pi(h_pi, basis, opt)(w), z0s, tol=pow(10.0, -10))

        ## see calc/stoh/l_3/res.d
        self.assertTrue(res.success)
        self.assertAlmostEqual(zs_opt[0], res.x[0])
        self.assertAlmostEqual(zs_opt[1], res.x[1])
        self.assertAlmostEqual(zs_opt[2], res.x[2])
        self.assertAlmostEqual((2.23413256581075-0.543249555891283j)*3, res.val)

    def test_opt_interface(self):
        opt_main(basis_type = 'STO',
                 basis_info = [(2, 0.6-0.4j, 'o')],
                 w0 = 1.0,
                 ws = np.linspace(0.5, 1.5, 11),
                 tol = pow(10.0, -5.0),
                 target = 'h_pi',
                 channel= '1s->kp',
                 dipole = 'length',
                 print_level = 1)
                 
                 
                 
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
