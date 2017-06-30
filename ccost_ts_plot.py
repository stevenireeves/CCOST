import matplotlib.pyplot as plt
import numpy as np
data = np.genfromtxt('CCOSTts.dat')
n = data.shape[0]
for i in range(1,n):
	plt.plot(data[0,:],data[i,:],label = "Oscillator" +str(i))
plt.xlabel('t')
plt.ylabel('Current')
plt.legend()
plt.show()
