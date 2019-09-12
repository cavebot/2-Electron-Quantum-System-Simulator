import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

data = np.loadtxt('pes_2D_C8.dat')

ecut = 0.6



#find number of independent energies along axis
counter1  = 1
k = 0 
while data[k,0] == data[k+1,0]:
	counter1 = counter1 + 1	
	k = k+1

counter2 = 0
#count energies below cutoff
for k in range(0, counter1):
	if data[k,1] < ecut:
		counter2 = counter2 + 1
	k = k+1


#Generate the photoelectron-spectrum for energies below a cutoff(e.g ecut=0.6)

cutdata = np.zeros((counter2**2,3))
#print(counter2)

length = 0 #length counter of the new population array after cutoff

k=0
with open("datafile.dat", "w") as output:

	for j in range(0, counter2):		
		for i in range(0,counter2):		
			if data[j*counter1 + i, 0] < ecut:
				if data[j*counter1 + i, 1] < ecut:
					cutdata[k,0] = data[j*counter1 + i, 0]
					cutdata[k,1] = data[j*counter1 + i, 1]
					cutdata[k,2] = data[j*counter1 + i, 2]
					output.write(str(cutdata[k,0])) 
					output.write("\t")  
					output.write(str(cutdata[k,1])) 
					output.write("\t")  
					output.write(str(cutdata[k,2]))
					output.write("\n")  
					k = k+1
output.close

cutdata2 = np.loadtxt('datafile.dat')
cutpop = cutdata2[:,2]

#generate population matrix for plotting as surface plot requires matrix input, not column
Z = np.zeros((counter2, counter2))
for i in range(0,counter2):
	for j in range(0, counter2):
		Z[i,j] = cutpop[i*counter2 + j]

#generate e1,e2 matrices as surface plot requires these to also be matrices :(

e1 = np.zeros((counter2, counter2))
e2 = np.zeros((counter2, counter2))

for j in range(0, counter2):
	for i in range(0, counter2):
		e1[j,i] = cutdata2[j*counter2 + i, 0] 
		e2[j,i] = cutdata2[j*counter2 + i, 1]


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel('$\epsilon_1 (s.a.u)$', fontsize=12, rotation=150)
ax.set_ylabel('$\epsilon_2 (s.a.u)$', fontsize=12, rotation=150)
ax.set_zlabel('population', fontsize=12, rotation=150)
ax.plot_surface(e1, e2, Z, cmap="autumn_r", lw=0.5, rstride=1, cstride=1, alpha=0.5)
ax.contour(e1, e2, Z, 10, lw=3, cmap="winter", linestyles="solid", offset=-0)
#ax.contour(e1, e2, Z, 10, lw=3, colors="k", linestyles="solid")
plt.show()



