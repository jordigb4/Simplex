from Simplex import Simplex
from Parser import parse

c, A, b  = parse(f'Datos/Datos23_4.txt')
sim = Simplex(c, A, b, debug_file= f'debug_Datos23_4')
sim.solve()

