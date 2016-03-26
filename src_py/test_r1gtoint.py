import numpy as np
from r1gtoint import *

import unittest
from minieigen import *

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
        res = self.gtos.calc_mat()
        self.assertAlmostEqual(1.0, res["s"][0, 0])
        self.assertAlmostEqual(1.0, res["s"][1, 1])

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

 
if __name__ == '__main__':
    unittest.main()
