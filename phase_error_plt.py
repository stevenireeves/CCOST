import matplotlib.pyplot as plt
import numpy as np
data = np.genfromtxt('ccost_phase_error.dat')
n = data.shape[0];
N = np.linspace(3,n+3,num = n);
ideal = data[0,0]*3/(N);
plt.plot(N,np.log10(data[:,0]), label = 'Average')
plt.plot(N,np.log10(data[:,1]), label = 'Min')
plt.plot(N,np.log10(data[:,2]), label = 'Max')
plt.plot(N,np.log10(ideal), label = 'Ideal')
plt.xlabel('N')
plt.ylabel('Phase Error')
plt.legend()
plt.show()
