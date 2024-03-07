import numpy as np
from Parser import parse


class Simplex:

    def __init__(self, c, A, b, __fase1 = False) -> None:
        self.__fase1 = __fase1
        self.__c = c
        self.__A = A
        self.__b = b

        self.__m = A.shape[0]; self.__n = A.shape[1]
        self.__dB = np.zeros(self.__m, dtype=int); self.__dN = np.zeros(self.__n - self.__m, dtype=int)

    def __isoptimal(self):
        pass
    
    def __get_direction(self):
        pass

    def __get_theta(self):
        pass

    def __update_values(self):
        pass

    def solve(self):

        #1. Find initial Basic Feasible Solution
        if not self.fase_1:
            self.__i_B, self.__i_N = self.__fase_1() 
        else:
            self.__i_N, self.__i_B = np.arange(1 , self.__n - self.__m + 1), np.arange(self.__n - self.__m, self._n + 1)

        B, c_B = A[self.__i_B], c[self.__i_B]
        A_N, c_N = A[self.__i_N], c[self.__i_N]
        
        #2. BFS optimal or choose entry variable

    def __fase_1(self):
        
        c = np.append(np.zeros(self.__n, dtype=int), np.ones(self.__m, dtype=int))

        identity_matrix = np.eye(self.__m)
        A = np.hstack((self.__A, identity_matrix))

        b = self.__b

        i_b, i_n = Simplex().solve(c, A, b, __fase1 = True)

        return i_b, i_n
        
c, A, b = parse('Datos/Datos5_2.txt')
Simplex = Simplex()
Simplex.solve(c, A, b)