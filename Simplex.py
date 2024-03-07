import numpy as np
from Parser import parse


class Simplex:

    def __init__(self, __fase1 = False) -> None:
        self.__fase1 = __fase1
    
    def __isoptimal(self):
        pass

    
    def __get_direction(self):
        pass

    def __get_theta(self):
        pass

    def __update_values(self):
        pass

    def solve(self, c, A, b):

        self.__c = c
        self.__A = A
        self.__b = b

        self.__m = A.shape[0]; self.__n = A.shape[1]
        self.__db = np.zeros(self.__m, dtype=int)
        self.__dn = np.zeros(self.__n - self.__m, dtype=int)

        if not self.fase_1:
            self.__i_b,self.__i_n =self.__fase_1()

        else:



            pass

    def __fase_1(self):
        
        c = np.append(np.zeros(self.__n, dtype=int), np.ones(self.__m, dtype=int))

        identity_matrix = np.eye(self.__m)
        A = np.hstack((self.__A, identity_matrix))

        b = self.__b

        i_b, i_n = Simplex().solve(c, A, b, __fase1 = True)

        return i_b,i_n
        
c, A, b = parse('Datos/Datos5_2.txt')
Simplex = Simplex()
Simplex.solve(c, A, b)