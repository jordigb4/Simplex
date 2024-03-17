import numpy as np
from Parser import parse


class Simplex:

    def __init__(self, c, A, b, phase1 = False, debug_file = 'debug'):
        """
        Initialize the data assocaited to an optimization problem
        """
        self._phase1 = phase1
        self._debug_file = debug_file
        self._iteration = 0

        self._m = A.shape[0]; self._n = A.shape[1]
        assert c.shape[0] == self._n and b.shape[0] == self._m, "Data dimensions are not correctly declared"
        self._c = c
        self._A = A
        self._b = b
        

    def _isoptimal(self):
        """
        Method to check if the basic feasible solution of a iteration is the optimal of the linear problem
        
        Pre-conditions: the working data is associated to a basic feasible solution

        Post-conditions: returns True if reduced costs are >= 0; otherwise fixes entry variable and returns False
        """

        self._r = self._c_N - (self._c_B @ (self._B_inv @ self._A_N))
        self._r = np.asarray(self._r).reshape(-1)

        if np.all(self._r >= 0):
            return True
        
        else: 
            negative_r = np.where(self._r < 0)[0]
            self._q = negative_r[0]  #considering that self.i_N is sorted (Blund)
            return False 


    def _get_direction(self):
        """
        Method to caluclate a basic feasible direction, associated to entry variable x_q

        Pre-conditions: local variable q can't be None or greater to (n-m)

        Post-conditions: declares an array with resulting DBF coordinates for basic variables;
                         if this array is >=0, the problem is not bounded and the program execution finishes
        """

        self._d_B = -self._B_inv @ np.asarray(self._A_N[:,self._q]).reshape(-1)
       
        if np.all(self._d_B >= 0):

            with open(f"{self._debug_file}.txt", "a") as debug_file:
                debug_file.write(f'The Problem is not bounded, it always can improve\n')
            raise Exception("The Problem is not bounded, it always can improve")
        

    def _get_theta(self):
        """
        Method to get the step length to reach the BFS that includes i_N[q] variable

        Pre-conditions: some element of d_B must be negative

        Post-conditions: declares the corresponding theta (step length) and exit variable associated
        """

        d_B = self._d_B
        theta = float('inf'); p = None
        for i, d_i in enumerate(d_B):   
            if d_i < 0:
                quotient =  -self._x_B[i] / d_i
                if quotient == theta: #Regla de Blund
                    p = p if self._i_B[p] < self._i_B[i] else i
                elif quotient < theta:
                    p = i
                    theta = quotient
                
        self._theta = theta
        self._p = p


    def _update_values(self):
        """
        Method to make the necessary updates of resulting BFS

        Pre-conditions: entry variable and DBF associated are calculated, also step length and its correponding exit variable
        
        Post-conditions: updates all variable related to a single BFS
        """

        self._B_inv = self._update_inverse()
        self._i_B[self._p], self._i_N[self._q] = self._i_N[self._q], self._i_B[self._p]
        self._i_N.sort()
        
        self._x_B += (self._theta * self._d_B)
        self._x_B[self._p] = self._theta
        
        self._A_N =  self._A[:,self._i_N]
        self._c_B, self._c_N = self._c[self._i_B], self._c[self._i_N]

        self._z += self._theta * self._r[self._q]


    def _update_inverse(self):
        """
        Method to update the basic's matrix inverse with reduced computational cost

        Pre-conditions: basic's matrix inverse, DBF and entry variable declarated must be of the same and previous simplex iteration

        Post-conditions: returns the inverse associated to the resulting BFS of the previous simplex iteration
        """
        
        E = np.eye(self._m)

        d_B = self._d_B
        d_Bp = d_B[self._p]
        n_p = (-d_B) / d_Bp
        n_p[self._p] = -1 / d_Bp

        E[:, self._p] = n_p
        return (E @ self._B_inv)
    
    
    def _clean_phase1(self):
        """
        Method to assert that all basic variables are from the orginal problem in phase 1
        """

        return np.all(self._i_B < self._n - self._m)
    
        
    def solve(self):
        """
        Method to solve the optimization problem associated to c, A and b

        Pre-conditions: problem data is already initialized and have the proper dimensions

        Post-conditions: if the problem corresponds to phase1 returns the BFS found and the basic matrix inverse associated,
                         so as to continue its updating; otherwise, returns the optimal solution, i.e. the basic variables and its value
        """

        #1. Find initial Basic Feasible Solution
        if not self._phase1:
            self._i_B, self._B_inv = self._phase_1()

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
                debug_file.write(f'Iteration {self._iteration}: p = {self._p}, q = {self._q}, B(p) = {self._i_B[self._p]}, theta* = {np.round(self._theta, 4)}, z = {np.round(self._z, 4)} \n')

        if self._phase1:
            if np.isclose(self._z, 0, atol=1e-10):
                while self._isoptimal() and not self._clean_phase1():
                    candidates = np.where(self._r == 0)[0]
                    self._q = candidates[0]
                    self._get_direction()
                    self._get_theta()
                    self._update_values()
                    self._iteration +=1

                with open(f"{self._debug_file}.txt", "a") as debug_file:
                    debug_file.write(f'Iteration {self._iteration}: p = {self._p}, q = {self._q}, B(p) = {self._i_B[self._p]}, theta* = {np.round(self._theta, 4)}, z = {np.round(self._z, 4)} \n')
                with open(f"{self._debug_file}.txt", "a") as debug_file:
                    debug_file.write(f'Basic feasible solution found, iteration {self._iteration} \n')

                return self._i_B, self._B_inv
            
            with open(f"{self._debug_file}.txt", "a") as debug_file:
                    debug_file.write('The problem is not feasible, artificial variables take value at optimal\n')
            raise Exception(f'The problem is not feasible')
        
        with open(f"{self._debug_file}.txt", "a") as debug_file:
                    debug_file.write(f'Optimal solution found, iteration {self._iteration}, z = {self._z} \n')
                    debug_file.write('End of primal simplex\n\n')
                    debug_file.write('Optimal solution: \n\n')
                    debug_file.write(f'vb = {self._i_B} \n')
                    debug_file.write(f'xb = {np.round(self._x_B, 3)} \n')
                    debug_file.write(f'z = {self._z} \n')
                    debug_file.write(f'r = {np.round(self._r, 2)}')

        return self._i_B, np.round(self._x_B, 3)
    

    def _phase_1(self):
        """
        Method to prepare data to phase 1 and solve the resulting linear problem
        """

        with open(f"{self._debug_file}.txt", 'w') as file:
            pass #cleaning the file
        with open(f"{self._debug_file}.txt", "a") as debug_file:
            debug_file.write("Start primal simplex with Bland's rule\n")
            debug_file.write('Phase I\n')

        c = np.append(np.zeros(self._n, dtype=int), np.ones(self._m, dtype=int))
        identity_matrix = np.eye(self._m)
        A = np.hstack((self._A, identity_matrix))
        b = self._b
        
        phase1 =  Simplex(c, A, b, phase1 = True, debug_file = self._debug_file)
        i_b, B_inv = phase1.solve()
        return i_b, B_inv
        
c, A, b = parse('Datos/Datos5_1.txt')
sim = Simplex(c, A, b, debug_file='Debug_Datos')
sim.solve()