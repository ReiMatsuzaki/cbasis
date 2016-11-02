import unittest
from r1basis import *
import scipy.special as sp
import numpy as np

class Test_driv_basis(unittest.TestCase):
    def setUp(self):
        pass

    def test_erfc(self):
        # in the application in STO-GTO integration,
        # the argument of erfcx becomes a/(2.0*sqrt(b)),
        # where a and b is orbital exponent of STO and GTO.
        # a ~ 1.0, 0.5 or 1/3
        # b ~ 0.01 ~ 30.0
        # So argument becomes
        # 5 ~ 0.027

        for theta in [0, 1, 45, 90, 270]:
            t = theta * np.pi / 180.0
            for a in [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]:
                z = a * np.exp(1.0j * a)
                self.assertAlmostEqual(sp.erfc(z),  erfc(z),  places=15)
                self.assertAlmostEqual(sp.erfcx(z), erfcx(z), places=15)

if __name__ == '__main__':
    unittest.main()
