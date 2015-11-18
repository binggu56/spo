import numpy as np 

import pylab as plt 

dat = np.genfromtxt(fname='fft.dat')

plt.plot(dat[:,0],dat[:,1],lw=2)

plt.show() 

 	
