from linspace import *
from set_l2func import *

import numpy as np

class BasisSet:
    def __init__(self, fs):
        self.fs = fs

    def matrix(self, op):
        return np.array([[cip(a, op(b)) for a in self.fs] for b in self.fs])

    def vector(self, b):
        return np.array([cip(a, b) for a in self.fs])

    def func(self, cs):
        return linear_combination(cs, self.fs)
