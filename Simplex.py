import numpy as np
from Parser import parse


class Simplex:

    def __init__(self, c, A, b, fase1 = False, debug_file = 'debug') -> None:    

        self.__fase1 = fase1
        self._m = A.shape[0]; self._n = A.shape[1]
        self._c = c
        self._A = A
        self._b = b
        self._debug_file = debug_file
        self._iteration = 0

    def _isoptimal(self):

        self._r = self._c_N - (self._c_B @ (self._B_inv @ self._A_N))
        self._r = np.asarray(self._r).reshape(-1)

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

        self._z += self._theta * self._r[self._q]


    def _update_inverse(self):
        
        E = np.eye(self._m)

        d_B = self._d_B
        d_Bp = d_B[self._p]
        n_p = (-d_B) / d_Bp
        n_p[self._p] = -1 / d_Bp

        E[:, self._p] = n_p
        return (E @ self._B_inv)
    
    def _clean_phase1(self):

        return np.all(self._i_B < self._n - self._m)
        
    def solve(self):

        #1. Find initial Basic Feasible Solution
        if not self.__fase1:
            self._i_B, self._B_inv = self._fase_1()
            with open(f"{self._debug_file}.txt", "a") as debug_file:
                debug_file.write('Fase II \n')
            self._B = self._A[:,self._i_B]
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

            self._iteration +=1
            with open(f"{self._debug_file}.txt", "a") as debug_file:
                debug_file.write(f'Iteració {self._iteration}: p = {self._p}, q = {self._q}, B(p) = {self._i_B[self._p]}, theta* = {np.round(self._theta, 4)}, z = {np.round(self._z, 4)} \n')

        if self.__fase1:

            if np.isclose(self._z, 0, atol=1e-10):
                while self._isoptimal() and not self._clean_phase1():
                    candidates = np.where(self._r == 0)[0]
                    self._q = candidates[0]
                    self._get_direction()
                    self._get_theta()
                    self._update_values()
                    self._iteration +=1
                with open(f"{self._debug_file}.txt", "a") as debug_file:
                debug_file.write(f'Iteració {self._iteration}: p = {self._p}, q = {self._q}, B(p) = {self._i_B[self._p]}, theta* = {np.round(self._theta, 4)}, z = {np.round(self._z, 4)} \n')
                with open(f"{self._debug_file}.txt", "a") as debug_file:
                    debug_file.write(f'Solució bàsica factible trobada, iteració {self._iteration} \n')
                return self._i_B, self._B_inv
            
            #artificial variables take value at optimal
            raise Exception("The Problem is not feasible")
        
        with open(f"{self._debug_file}.txt", "a") as debug_file:
                    debug_file.write(f'Solució òptima trobada, iteració {self._iteration}, z = {self._z} \n')
                    debug_file.write('Fi simplex primal \n\n')
                    debug_file.write('Solució òptima: \n\n')
                    debug_file.write(f'vb = {self._i_B} \n')
                    debug_file.write(f'xb = {np.round(self._x_B, 3)} \n')
                    debug_file.write(f'z = {self._z} \n')
                    debug_file.write(f'r = {np.round(self._r, 2)}')

    def _fase_1(self):

        with open(f"{self._debug_file}.txt", 'w') as file:
            pass #cleaning the file
        with open(f"{self._debug_file}.txt", "a") as debug_file:

            debug_file.write('Inici simplex primal amb regla de Bland \n')
            debug_file.write('Fase I\n')

        c = np.append(np.zeros(self._n, dtype=int), np.ones(self._m, dtype=int))
        identity_matrix = np.eye(self._m)
        A = np.hstack((self._A, identity_matrix))
        b = self._b
        
        fase1 =  Simplex(c, A, b, fase1 = True, debug_file = self._debug_file)
        i_b, B_inv = fase1.solve()
        return i_b, B_inv
        
c, A, b = parse('Datos/Datos5_1.txt')
sim = Simplex(c, A, b, debug_file='Debug_Datos_5_1')
sim.solve()