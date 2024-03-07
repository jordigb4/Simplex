import numpy as np
from Parser import parse


class Simplex:

    def __init__(self, c, A, b) -> None:
    
        self._m = A.shape[0]; self._n = A.shape[1]
        self._c = np.insert(c, 0, None)
        self._A = np.c_[np.full(self.__m, None), A]
        self._b = b

        self._dB = np.zeros(self._m, dtype=int); self._dN = np.zeros(self._n - self._m, dtype=int)

    def _isoptimal(self):
        
        self._r_N = self._c_N - self._c_B
    
    def _get_direction(self):
        pass

    def _get_theta(self):
        pass

    def _update_values(self):
        pass

    def _update_inverse(self):
        pass

    def solve(self, __fase1 = False):

        #1. Find initial Basic Feasible Solution
        if not __fase1:
            self._i_B, self._i_N = self.__fase_1() 
        else:
            self._i_N, self._i_B = np.arange(1 , self._n - self._m + 1), np.arange(self._n - self._m, self._n + 1)

        self._B, self._c_B = A[self._i_B], c[self._i_B]
        self._A_N, self._c_N = A[self._i_N], c[self._i_N]
        
        #2. BFS optimal or choose entry variable

    def _fase_1(self):
        
        c = np.append(np.zeros(self._n, dtype=int), np.ones(self._m, dtype=int))

        identity_matrix = np.eye(self._m)
        A = np.hstack((self._A, identity_matrix))

        b = self._b

        i_b, i_n = Simplex().solve(c, A, b, __fase1 = True)

        return i_b, i_n
        
c, A, b = parse('Datos/Datos5_2.txt')
Simplex = Simplex()
Simplex.solve(c, A, b)