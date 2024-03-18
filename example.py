from Simplex import Simplex
from Parser import parse

c, A, b  = parse(f'Datos/Datos5_2.txt')
sim = Simplex(c, A, b, debug_file= f'debug_Datos5_2')
sim.solve()

