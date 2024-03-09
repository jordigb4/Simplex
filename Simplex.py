import numpy as np
from Parser import parse


class Simplex:

    def __init__(self, c, A, b, fase1 = False) -> None:

        self.__fase1 = fase1
        self._m = A.shape[0]; self._n = A.shape[1]
        self._c = c
        self._A = A
        self._b = b
        self.__iteration =0

    def _isoptimal(self):
        
        print(f"Simplex {'fase 1' if self.__fase1 else 'primal' }. Iteració: {self.__iteration}\n\n")
        print(f"x_B = {self._x_B}")
        print(f"i_B = {self._i_B}")
        print(f"B_inv = {self._B_inv}")
        print(f"z = {self._z}")
        print("====================================================\n\n")
        self.__iteration +=1

        self._r = self._c_N - (self._c_B @ (self._B_inv @ self._A_N))
        
        if np.all(self._r >= 0):
            return True
        
        else:
            negative_r = np.where(self._r < 0)[0]
            self._q = negative_r[0]  #considering that self.i_N is sorted (Blund)
            return False 


    def _get_direction(self):

        self._d_B = -self._B_inv @ np.asarray(self._A_N[:,self._q]).reshape(-1)
       
        if np.all(self._d_B >= 0):
            raise Exception("The Problem is not bounded, it always can improve")
        
    def _get_theta(self):
        
        d_B = self._d_B
        theta = float('inf'); p = None
        for i, d_i in enumerate(d_B):   
            if d_i < 0:
                quocient =  -self._x_B[i] / d_i
                if quocient == theta: #Regla de Blund
                    p = p if self._i_B[p] < self._i_B[i] else i
                elif quocient < theta:
                    p = i
                    theta = quocient
                
        self._theta = theta
        self._p = p

    def _update_values(self):

        self._B_inv = self._update_inverse()
        self._i_B[self._p], self._i_N[self._q] = self._i_N[self._q], self._i_B[self._p]
        self._i_N.sort()
        
        self._x_B += (self._theta * self._d_B)
        self._x_B[self._p] = self._theta
        
        self._A_N =  self._A[:,self._i_N]
        self._c_B, self._c_N = self._c[self._i_B], self._c[self._i_N]

        self._z += self._theta * self._r[ :, self._q]

    def _update_inverse(self):
        
        E = np.eye(self._m)

        d_B = self._d_B
        d_Bp = d_B[self._p]
        n_p = (-d_B) / d_Bp
        n_p[self._p] = -1 / d_Bp

        E[:, self._p] = n_p
        return (E @ self._B_inv)
        
    def solve(self):

        #1. Find initial Basic Feasible Solution
        if not self.__fase1:
            self._i_B, self._B_inv = self._fase_1()
            self._B = self._A[self._i_B]
            self._i_N = sorted(np.setdiff1d(np.arange(self._n), self._i_B))
        else:
            self._i_N, self._i_B = np.arange(0 , self._n - self._m), np.arange(self._n - self._m, self._n)
            self._B = self._B_inv = np.eye(self._m)

        
            
        self._c_B, self._c_N = self._c[self._i_B], self._c[self._i_N]
        self._A_N  = self._A[:,self._i_N]
        self._x_B = self._B_inv @ b; self._x_N = np.zeros(self._n - self._m, dtype=int)

        self._z = self._c_B @ self._x_B
            
            
        #2. BFS optimal or choose entry variable q and continue  
        while not self._isoptimal():
            #3. Get BFD associated to q
            self._get_direction()
            #4. Get maximum pass length and exit variable p
            self._get_theta()
            #5. Update structures
            self._update_values()
            
             
        if self.__fase1:
            if self._z == 0:
                if np.all(self._i_B < self._m):
                    return self._i_B, self._B_inv
                else:
                    #Degenerated solution of Problem Phase1
                    i_artificial = np.where(self._i_B >= self._m)[0]
                    i_original = np.where(self._i_N < self._m)[0]
                    for i in i_artificial:
                        j = i_original.pop()
                        self._i_B[i] = self._i_N[j]
                    return self._i_B, self._B_inv
            #artificial variables take value at optimal
            raise Exception("The Problem is not feasible")
                     
    def _fase_1(self):
        
        c = np.append(np.zeros(self._n, dtype=int), np.ones(self._m, dtype=int))
        identity_matrix = np.eye(self._m)
        A = np.hstack((self._A, identity_matrix))
        b = self._b
        
        fase1 =  Simplex(c, A, b, fase1 = True)
        i_b, B_inv = fase1.solve()
        return i_b, B_inv
        
c, A, b = parse('Datos/toy.txt')
sim = Simplex(c, A, b)
sim.solve()