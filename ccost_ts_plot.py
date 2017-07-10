import matplotlib.pyplot as plt
import numpy as np
data = np.genfromtxt('CCOSTts.dat')
u = np.zeros(data.shape[1])
n = data.shape[0]
f, axarr = plt.subplots(2, sharex =True)
for i in range(1,n):
	axarr[0].plot(data[0,:],data[i,:],label = "Oscillator" +str(i))
	u += data[i,:]
axarr[1].plot(data[0,:],u)
plt.xlabel('t')
#axarr[0].ylabel('Current')
#axarr[1].ylabel('Summed Signal')
plt.legend()
plt.show()
