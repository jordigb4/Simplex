import numpy as np

def parse(path):
    """
    Function to read optimization problem data from a file

    Pre-conditions: path is a valid path to a readable file. Problem's parameters are written in the following lines
                    of pattern '<parameter>=', having parameter = {c, A, b}. Only numbers must be written.

    Post-condition: creates from data in path a vector with cost coefficients c, antoher with right-hand side of 
                    constraints b, and a matrix with coefficients of variables for each constraint A
    """
    
    with open(path, 'r') as file:

        lines = [line.strip() for line in file.readlines() if line.strip()]

    # Dict to store the parsed data
    variables = {}

    current_variable = None
    for line in lines:

        if line.endswith('='):
            # Get variable name and add to dict
            current_variable = line[:-1]  
            variables[current_variable] = []    

        elif current_variable is not None:
            # Extract and store in variable key dict
            values = [int(val) for val in line.split()]
            if current_variable == 'A':
                variables[current_variable].append(values)
            else:
                variables[current_variable].extend(values)

    return np.array(variables['c']), np.matrix(variables['A']), np.array(variables['b'])
