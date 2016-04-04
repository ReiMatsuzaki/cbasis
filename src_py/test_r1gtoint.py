import numpy as np
#from r1gtoint import *
from r1gtoint import *

import unittest
from minieigen import *

class Test_Eigen(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_print_matrixxi(self):
        ivec = VectorXi.Zero(3)
        ivec[0] = 1; ivec[1] = 2; ivec[2] = 3;
        print_vectorxi(ivec)
        print_vectorxi([4, 3, 3])

class Test_r1gtos(unittest.TestCase):
    def setUp(self):
        self.gtos = R1GTOs(0)
        self.gtos.add( 2, 1.1)
        self.gtos.add(3, [1.2, 1.3-0.1j])

    def test_r1sto_vec(self):
        stos = vector_R1STO()
        stos.append(R1STO(1.1, 2, 1.2))
        stos.append(R1STO(1.1-0.2j, 3, 1.3))
        self.assertAlmostEqual(1.1-0.2j, stos[1].c)

    def test_add(self):
        self.assertEqual(2, self.gtos.basis(0).n)
        self.assertAlmostEqual(1.1, self.gtos.basis(0).z)
        self.assertAlmostEqual(1.2, self.gtos.basis(1).z)
        self.assertAlmostEqual(1.3-0.1j, self.gtos.basis(2).z)

    def test_matrix(self):
        self.gtos.calc_mat()
        self.assertAlmostEqual(1.0, self.gtos.mat("s")[0, 0])
        self.assertAlmostEqual(1.0, self.gtos.mat("s")[1, 1])

    def test_vector(self):
        stos = R1STOs()
        stos.add(1.1, 2, 1.2)
        stos.add(1.1-0.2j, 3, 1.3)
        res = self.gtos.calc_vec(stos)

    def test_print(self):
        s = R1STO(1.0, 2, 1.1)
        print s

    def test_at_r(self):
        gs = R1GTOs(0)
        gs.add(1.1, 2, 1.2)
        ys = gs.at_r([1.3], [1.1])[0]
        print ys

    def test_opt_alpha_shift(self):
        gs = R1GTOs(1)
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
        
        driv = R1STOs(); driv.add(2.0, 2, 1.0)
        
        #opt_idx = np.array([7.0, 8, 9, 10, 11, 12, 13, 14]);
        #opt_idx = np.array([7, 8, 9, 10, 11, 12, 13, 14]);
        opt_idx = [7, 8, 9, 10, 11, 12, 13, 14];

        h= 0.0001;
        eps = 0.000001;
        z_shift = 0.0-0.02j;
        convq = False;
        alpha = 0.0;
        (convq, alpha, z_shift) = opt_alpha_shift(driv, opt_idx,
                                                  0.5,
		                                  h, 100, eps, 10.0**(-10), 
		                                  gs, z_shift)
        self.assertTrue(convq)
        shift_ref = -0.00293368-0.0204361j
        self.assertAlmostEqual(shift_ref, z_shift, places=4);
        alpha_ref = -5.6568937518988989+1.0882823480377297j
        self.assertAlmostEqual(alpha_ref, alpha)

    def test_solve_alpha(self):
        gs = R1GTOs(1)
        zs = [0.463925,
              1.202518,
              3.379649,
              10.6072,
              38.65163,
              173.5822,
              1170.498,
              0.16934112166516593 - 0.01j,
              0.08989389391311804 - 0.01j,
              0.05561087391349172 - 0.01j,
              0.03776599632952126 - 0.01j,
              0.02731159914174668 - 0.01j,
              0.02066585522406014 - 0.01j,
              0.01618060242100465 - 0.01j,
              0.01301156966796773 - 0.01j]
        gs.add(2, zs)
        driv = R1STOs(); driv.add(2.0, 2, 1.0)
        ref = calc_alpha(driv, gs, 0.57, 10.0**(-7))
        cs = solve_alpha(driv, gs, 0.57, 10.0**(-7))
        gs.calc_vec(driv)
        calc= sum([c*m for (c, m) in zip(np.array(cs), gs.vec("m"))])
        self.assertAlmostEqual(ref, calc)
        
 
if __name__ == '__main__':
    unittest.main()
