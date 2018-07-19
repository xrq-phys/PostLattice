#!/usr/bin/env python

import numpy as np

class double_occ:
    '''Double Occupacy'''

    def __init__(self, system):
        self.value = float(0.)
        self.system = system

    def measure(self, i, si, j, sj, k, sk, l, sl, x):
        if (i == j and i == k and i == l and si != sk and sj != sl and si == sj): 
            self.value += x / 2.

class spin_struct:
    '''Spin Structure'''

    def __init__(self, system, points):
        self.system = system
        self.points = points
        self.values = np.zeros(np.size(points), complex)
        self.paulix = np.array([ [ 0,  1 ], [ 1,  0 ] ])
        self.pauliy = np.array([ [ 0, -1 ], [ 1,  0 ] ])  # i is expressed as i^2=-1 hence S^2 = Sx2 + Sz2 - Sy2
        self.pauliz = np.array([ [ 1,  0 ], [ 0, -1 ] ])

    def measure(self, i, si, j, sj, k, sk, l, sl, x):
        if (i == j and k == l):
            spin_part = self.paulix[si][sj] * self.paulix[sk][sl] \
                      + self.pauliz[si][sj] * self.pauliz[sk][sl] \
                      - self.pauliy[si][sj] * self.pauliy[sk][sl]

            for ip in range(0, np.size(self.points, 0)):
                space_part = 1 / (3. * self.system.n ** 2) \
                           * np.exp(1j * np.inner(self.points[ip], self.system.r(i) - self.system.r(k)))
                self.values[ip] += space_part * spin_part * x

class sc_corr:
    '''Superconducting Correlation'''

    def __init__(self, system, rc_list):
        self.system = system
        self.rc_list = np.array(rc_list, int)
        self.values = np.zeros(np.size(self.rc_list), float)

    def measure(self, nu, snu, mu, smu, xi, sxi, eta, seta, x):
        sign = 1
        r_d = self.system.r(mu) - self.system.r(xi)
        r_to = self.system.r(nu) - self.system.r(mu)
        r_from = self.system.r(xi) - self.system.r(eta)

        r_m = 10000
        for i in [ [0, 0], [0, 1], [0, -1], [1, 0], [1, 1], [1, -1], [-1, 0], [-1, 1], [-1, -1] ]:
            r_d = self.system.r(mu) - self.system.r(xi) + np.array(i) * self.system.a
            if np.dot(r_d, r_d) < r_m:
                r_m = np.dot(r_d, r_d)

        if (snu != smu and sxi != seta and self.system.nn[nu][mu] and self.system.nn[xi][eta]):

            if (smu != sxi):
                sign = -sign
            if (np.dot(r_from, r_to) == 0): 
                sign = -sign
            
            for i in range(np.size(self.rc_list)):
                if (r_m == self.rc_list[i]):
                    self.values[i] += x * sign / (2. * self.system.n)
