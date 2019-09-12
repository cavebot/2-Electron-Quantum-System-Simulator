from scipy.integrate import simps
import numpy as np

file = "pes_2D_C5.dat"

print(file)
x = np.loadtxt(file, usecols=1)
y = np.loadtxt(file, usecols=1)
z = np.loadtxt(file, usecols=2)

grid = 58
x = x[0:grid]
y = y[0:grid]

R  = np.zeros((grid, grid))
for j in range(0, grid-1):
	for i in range(0, grid-1):
		R[i,j] = z[i + j*grid]
		
print(len(R))

b = simps(simps(R, x), y)
print(b)

#0.85134099743259539
#import sympy
#xx, yy = sympy.symbols('x y')
#sympy.integrate(sympy.cos(xx)**4 + sympy.sin(yy)**2, (xx, 0, 1), (yy, 0, 1)).evalf()
#0.851349922021627
