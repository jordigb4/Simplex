import numpy as np

def parse(path):
    with open(path, 'r') as file:

        lines = [line.strip() for line in file.readlines() if line.strip()]

    # Variables to store the parsed data
    variables = {}

    # Read line by line
    current_variable = None
    for line in lines:
        if line.endswith('='):
            current_variable = line[:-1]  # What is the variable name?
            variables[current_variable] = []    
        elif current_variable is not None:
            # Extract and store in variable
            values = [int(val) for val in line.split()]
            if current_variable == 'A':
                variables[current_variable].append(values)
            else:
                variables[current_variable].extend(values)


    return np.array(variables['c']), np.matrix(variables['A']), np.array(variables['b'])

c, A, b = parse('Datos/Datos5_2.txt')

'''theta = float('inf')
        for i in negative_directions.sort():
            quocient =  -self._x_B[i] / d_B[i]
            if quocient < theta:
                theta = quocient
                p = i'''

'''candidates = -self._x_B[negative_directions] / d_B[negative_directions]
        min_value = np.min(candidates)
        if np.count_nonzero(candidates == min_value)>=2:
            self._theta = np.min(candidates)
            self._p = np.argmin(self._i_B[negative_directions])
    
        else:

            self._theta, self._p = np.min(candidates), negative_directions[np.argmin(candidates)]'''


