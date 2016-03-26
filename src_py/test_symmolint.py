import numpy as np
from symmolint import *

import unittest
import minieigen as me
import scipy.linalg as la

class Test_minieigen(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_print(self):
        v = me.VectorXc([5, 6])
        x = me.MatrixXc([[1, 2], [3, 4]])
        xv = x*v

class Test_ceig(unittest.TestCase):
    def setUp(self):
        pass

    def test_ceig(self):
        h = me.MatrixXc(
            [[1.0+0.1j,  0.2,   0.1-0.1j ],
             [0.2,       0.3,   0.4      ],
             [0.1-0.1j,  0.4,   1.0-0.2j ]])
        s = me.MatrixXc(
            [[1.0,       0.2,   0.2+0.00j ],
             [0.2,       1.0,   0.1       ],
             [0.2+0.00j, 0.1,   1.0       ]])
        (eigs, eigvec) = ceig(h, s)
        for i in range(h.cols()):
            lexp = h*eigvec.col(i)
            rexp = s*eigvec.col(i)*eigs[i]
            self.assertAlmostEqual(0.0, abs(lexp-rexp))

        for i in range(eigs.rows()-1):
            self.assertTrue(eigs[i].real < eigs[i+1].real)

            
class Test_BMatSet(unittest.TestCase):
    def setUp(self):
        pass

    def test_getset(self):
        cs = Cs()
        bmat = BMatSet(cs.order)
        s00 = me.MatrixXc.Zero(4, 4)
        s00[0, 1] = 1.1; s00[1, 0] = 1.2
        bmat.set_matrix("s", Cs_Ap(), Cs_Ap(), s00)
        s00_get1 = bmat["s", Cs_Ap(), Cs_Ap()]
        s00_get2 = bmat["s", Cs_Ap(), Cs_Ap()]
        s00_get1[0, 1] = 2.1
        self.assertAlmostEqual(1.2, s00_get2[1, 0])
        self.assertAlmostEqual(2.1, s00_get2[0, 1])
        
        
class Test_SymMolInt(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_symmetry_group(self):
        sym = Cs()
        self.assertEqual(2, sym.order)
        self.assertEqual("Cs", sym.name)

    def test_reduction_sets(self):
        """
        coef_iat_ipn = me.MatrixXc.Zero(2, 3)
        coef_iat_ipn[0, 1] = 1.1
        coef_iat_ipn[1, 0] = 1.2
        """
        coef_iat_ipn = [[0.0, 1.1], [1.2, 0.0]]
        rds1 = Reduction(Cs_Ap(), coef_iat_ipn)
        self.assertEqual(Cs_Ap(), rds1.irrep)
        self.assertAlmostEqual(1.1, rds1.coef_iat_ipn()[0, 1])
        self.assertAlmostEqual(1.2, rds1.coef_iat_ipn()[1, 0])

    def test_sub(self):
        sub = sub_two_sgto(Cs(), Cs_App(), [1,2,3], [2,4,5])
        self.assertAlmostEqual(1.0, sub.get_xyz_iat()[0, 0])
        self.assertAlmostEqual(1.0, sub.get_xyz_iat()[0, 1])
        self.assertAlmostEqual(2.0, sub.get_xyz_iat()[1, 0])
        self.assertAlmostEqual(2.0, sub.get_xyz_iat()[1, 1])
        self.assertAlmostEqual(3.0, sub.get_xyz_iat()[2, 0])
        self.assertAlmostEqual(-3.0, sub.get_xyz_iat()[2, 1])

        self.assertAlmostEqual(2.0, sub.get_zeta_iz()[0])
        self.assertAlmostEqual(4.0, sub.get_zeta_iz()[1])
        self.assertAlmostEqual(5.0, sub.get_zeta_iz()[2])

    def test_sub2(self):
        sub = (SubSymGTOs()
               .xyz((0, 0.1, 0.2))
               .ns((0, 1, 2))
               .xyz((2, 1, 0))
               .rds(Reduction(Cs_Ap(),  [[+1, +1]]))
               .rds(Reduction(Cs_App(), [[+1, -1]]))
               .zeta([1.1, 1.4, 2.1]))

        self.assertAlmostEqual(0.1, sub.get_xyz_iat()[1, 0])
        self.assertAlmostEqual(1.0, sub.get_xyz_iat()[1, 1])
        
    def test_s_center(self):
        gtos = SymGTOs(Cs())
        zs = [2.0**n-0.1j for n in range(-2, 2)]
        gtos.sub(sub_s(Cs_Ap(), (0, 0, 0), zs))
        gtos.atom([0, 0, 0], 1.1)
        gtos.setup()
        mat_set = gtos.calc_mat()

        cs = [1.0, 1.1, 1.2, 1.3]
        rs = [1.0, 2.0]
        ys = gtos.at_r_ylm(0, 0, Cs_Ap(), cs, rs)


class Test_H_atom(unittest.TestCase):
    def setUp(self):
        xatom = (0, 0, 0)
        Ap = Cs_Ap()
        App= Cs_App()
        self.gtos = (SymGTOs(Cs())
                     .sub(SubSymGTOs()
                          .xyz(xatom)
                          .ns((0, 0, 0))
                          .ns((0, 0, 1))
                          .rds(Reduction(Ap,  [[1, 0]]))
                          .rds(Reduction(App, [[0, 1]]))
                          .zeta([2.0**n for n in range(-10, 10)]))
                     .atom(xatom, 1.0))
        mat = self.gtos.calc_mat()
        
        s = mat["s", Ap, Ap]
        t = mat["t", Ap, Ap]
        v = mat["v", Ap, Ap]
        h = t + v
        (self.eigs, self.eigvecs) = ceig(h, s)
     
    def test_energy(self):
        for i in range(0, 4):
            n = i + 1
            ene = -0.5/(n*n);
            print ene, self.eigs[i]
            self.assertTrue(abs(ene-self.eigs[i])   < 0.0001);

    def test_wavefunc(self):
        c = self.eigvecs.col(0)
        c = self.gtos.correct_sign(0, 0, 0, c)
        rs = [1.1]
        (ys_calc, ds_calc) = self.gtos.at_r_ylm(0, 0, Cs_Ap(), c, rs)
        ys_refs = [2.0*r*np.exp(-r) for r in rs]
        self.assertTrue(abs(ys_calc[0] - ys_refs[0]) < 0.0001)

    def test_2p(self):
        App= Cs_App()
        mat = self.gtos.calc_mat()
        s = mat["s", App, App]
        t = mat["t", App, App]
        v = mat["v", App, App]        
        (eigs, eigenvecs) = ceig(t+v, s)

        self.assertTrue(abs(eigs[0]+0.125) < 0.000001)

        c = self.gtos.correct_sign(1, 0, App, eigenvecs.col(0))
        rs = [1.1]
        (ys_calc, ds_calc) = self.gtos.at_r_ylm(1, 0, App, c, rs)

        ys_refs = [(1.0/(2.0*np.sqrt(6.0)))*r*r*np.exp(-0.5*r) for r in rs]
        self.assertTrue(abs(ys_calc[0] - ys_refs[0]) < 0.0001)        
        

class Test_H2_plus(unittest.TestCase):
    def setUp(self):
        self.gtos = (SymGTOs(Cs())
                     .sub(SubSymGTOs()
                          .xyz((0, 0, 0.7))
                          .xyz((0, 0,-0.7))
                          .ns((0, 0, 0))
                          .ns((0, 0, 1))
                          .rds(Reduction(Cs_Ap(), [[1,0],
                                                   [1,0]]))
                          .rds(Reduction(Cs_Ap(), [[0,1],
                                                   [0,-1]]))
                          .zeta([2.0**n for n in range(-10,10)]))
                     .atom((0, 0, 0.7), 1.0)
                     .atom((0, 0,-0.7), 1.0))

    def test_energy(self):
        """
        reference energy is from 
        calc/ccolumbus/h2/theta_traj/lie_2s3p_4s6dcstong/out/theta10_10g/out
        in rcclsc
        """
        self.gtos.setup()
        mat = self.gtos.calc_mat()
        Ap = Cs_Ap()
        s = mat["s", Ap, Ap]
        t = mat["t", Ap, Ap]
        v = mat["v", Ap, Ap]
        h = t + v
        (eigs, eigvecs) = ceig(h, s)
        self.assertAlmostEqual(-1.284146, eigs[0], places=4,
                               msg = "eigs[0] = {0}".format(eigs[0]))


class Test_H_photoionization(unittest.TestCase):
    def calc_full(self):
        xatom = (0, 0, 0)
        Ap = Cs_Ap()
        App= Cs_App()
        zeta0 = [2.0**n for n in range(-10, 10)]
        zeta1 = [2.0**n-0.02j for n in range(-15, 5)]
        self.gtos = (SymGTOs(Cs())
                     .sub(SubSymGTOs()
                          .xyz(xatom)
                          .ns((0, 0, 0))
                          .rds(Reduction(Ap,  [[1]]))
                          .zeta(zeta0))
                     .sub(SubSymGTOs()
                          .xyz(xatom)
                          .ns((0, 0, 1))
                          .rds(Reduction(App, [[1]]))
                          .zeta(zeta1))
                     .atom(xatom, 1.0))
        self.gtos.setup()
        mat = self.gtos.calc_mat()
        h0 = mat["t", Ap, Ap] + mat["v", Ap, Ap]
        s0 = mat["s", Ap, Ap]
        (eigs0, eigvecs0) = ceig(h0, s0)
        ene0 = eigs0[0]
        c0 = self.gtos.correct_sign(0, 0, Ap, eigvecs0.col(0))
        self.assertAlmostEqual(ene0, -0.5, places=4)
        
        h1 = mat["t", App, App] + mat["v", App, App]
        s1 = mat["s", App, App]

        z10 = mat["z", App, Ap]
        w = 1.0
        m1 = z10 * c0
        c1 = la.solve(h1 - (w+ene0)*s1, m1)
        self.alpha_full = (m1 * c1).sum()

    def calc_part(self):
        xatom = (0, 0, 0)
        Ap = Cs_Ap()
        App= Cs_App()
        zeta0 = [2.0**n for n in range(-10, 10)]
        zeta1 = [2.0**n-0.02j for n in range(-15, 5)]
        self.gtos0 = (SymGTOs(Cs())
                      .sub(sub_s(Ap, xatom, zeta0))
                      .atom(xatom, 1.0))
        self.gtos1 = (SymGTOs(Cs())
                      .sub(sub_pz(App, xatom, zeta1))
                      .atom(xatom, 1.0))
        mat0 = self.gtos0.calc_mat()
        mat1 = self.gtos1.calc_mat()
        mat10= self.gtos1.calc_mat(self.gtos0, False)
        h0 = mat0["t", Ap, Ap] + mat0["v", Ap, Ap]
        s0 = mat0["s", Ap, Ap]
        (eigs0, eigvecs) = ceig(h0, s0)
        ene0 = eigs0[0]
        c0 = eigvecs.col(0)
        h1 = mat1["t", App, App] + mat1["v", App, App]
        s1 = mat1["s", App, App]
        z10 = mat10["z", App, Ap]
        m1 = z10 * c0
        w = 1.0
        c1 = la.solve(h1 - (w+ene0)*s1, m1)
        self.alpha_part = (m1 * c1).sum()
        
    def setUp(self):
        self.calc_full()
        self.calc_part()

    def test_1skp(self):
        # see calc/stoh/1skp/l_5 in rcclsc (E=1)
        alpha_ref = -1.88562800720386+0.362705406693342j
        self.assertAlmostEqual(self.alpha_full, alpha_ref, places=3)
        self.assertAlmostEqual(self.alpha_part, self.alpha_part)

        
if __name__ == '__main__':
    unittest.main()
        
