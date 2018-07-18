#!/usr/bin/env python

import numpy as np

class square_2d:
    '''2D Square Lattice'''

    def __init__(self, a):
        self.a  = np.array(a, int)
        self.n  = np.product(a)
        self.nn = np.zeros([ self.n, self.n ], int)

        for yi in range(0, a[0]):
            for xi in range(0, a[1]):
                self.nn[self.idx_2d(xi, yi)][self.idx_2d(xi, (yi + 1)        % a[0])] = 1
                self.nn[self.idx_2d(xi, yi)][self.idx_2d(xi, (yi - 1 + a[0]) % a[0])] = 1
                self.nn[self.idx_2d(xi, yi)][self.idx_2d((xi + 1)        % a[0], yi)] = 1
                self.nn[self.idx_2d(xi, yi)][self.idx_2d((xi - 1 + a[0]) % a[0], yi)] = 1

    def idx_2d(self, x, y):
        return y * self.a[1] + x

    def r(self, i):
        return np.array([ i % self.a[1], i / self.a[1] ], int)

