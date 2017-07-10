import matplotlib.pyplot as plt
import numpy as np
data = np.genfromtxt('ccost_phase_error.dat')
n = data.shape[0];
num = 3
for i in range(1,n):
	num += 2
N = np.linspace(3,num, n)
ideal = data[0,0]*3/(N)
plt.plot(N,np.log10(data[:,0]), label = 'Average')
plt.plot(N,np.log10(data[:,1]), label = 'Min')
plt.plot(N,np.log10(data[:,2]), label = 'Max')
plt.plot(N,np.log10(ideal), label = '1/N')
plt.xlabel('N')
plt.ylabel('log Phase Error')
plt.legend()
plt.show()
