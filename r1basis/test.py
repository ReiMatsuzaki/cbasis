import numpy as np
from r1basis import *

import unittest


"""
class Test_Eigen(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_print_matrixxi(self):
        ivec = VectorXi.Zero(3)
        ivec[0] = 1; ivec[1] = 2; ivec[2] = 3;
        print_vectorxi(ivec)
        print_vectorxi([4, 3, 3])

"""

class Test_first(unittest.TestCase):
    def setUp(self):
        pass

    def test_add(self):
        print add(1,2)


class Test_r1_linear_comb(unittest.TestCase):
    def setUp(self):
        pass

    def test_size(self):

        for (m, fs) in zip([1], [LC_STOs()]): 
            fs.add(1.1, 2, 1.3-0.2j)
            fs.add(0.1, 1, 1.2-0.2j)
            fs.add(0.2, 3, 1.4-0.1j)
            self.assertEqual(3, fs.size())
        
        
    def test_at_r(self):
        
        r = 1.4
        c0 = 1.1; n0 = 1; z0 = 0.35; 
        c1 = 1.2; n1 = 2; z1 = 0.3; 

        #for (m, fs) in zip([1,2], [LC_STOs(), LC_GTOs()]): 
        for (m, fs) in zip([1], [LC_STOs()]): 
            fs.add(c0, n0, z0)
            fs.add(c1, n1, z1)
            print fs.str()
            ys = fs.at_r([r])
            y_ref = (c0*r**n0*np.exp(-z0*r**m) + 
                     c1*r**n1*np.exp(-z1*r**m))
            self.assertAlmostEqual(y_ref, ys[0])

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
