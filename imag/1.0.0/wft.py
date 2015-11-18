import numpy as np 

import pylab as plt 

dat0 = np.genfromtxt(fname='wf0.dat')
dat1 = np.genfromtxt(fname='wft.dat')

#dat2 = np.genfromtxt(fname='/home/bing/dvr/1.0.2/wf_dat')
#plt.plot(dat0[:,0],dat0[:,1],lw=2,label='$|\psi(x,0)|^2$')
plt.plot(dat1[:,0],dat1[:,1],lw=2,label='$|\psi(x,t)|^2$')

#plt.plot(dat2[:,0],dat2[:,1]**2,lw=2,label='$|\psi(x,t)|^2$')

plt.legend()

plt.xlim(0,3)

plt.show() 

 	
