##!/usr/bin/python

import numpy as np
import pylab as plt
import cmath 
#import seaborn as sns 

#sns.set_context('poster')

plt.subplot(1,1,1)
t,x,y = np.genfromtxt(fname='corr',usecols=(0,1,2),unpack=True) 

#for x in range(1,data.shape[-1]):
    
#plt.plot(data[:,0],data[:,1],lw=2,label='Re')
#plt.plot(data[:,0],data[:,2],lw=2,label='Im')
#plt.plot(data[:,0],data[:,3],lw=2,label='$|C(t)|$')

#dat = np.genfromtxt(fname='../spo/1.0.2/corr')
#plt.plot(dat[:,0],dat[:,1],'--',label='Re, QM')
#plt.plot(dat[:,0],dat[:,2],'--',label='Im, QM')

#x = np.linspace(0,4,101) 
#y = []
#for i in range(0,101):
#!    y[i] = 
#for i in range(0,101):
#y = 0.3*np.exp(2j*np.pi*x)+0.5*np.exp(1.2j*2.0*np.pi*x)+0.2*np.exp(6j*2.0*np.pi*x)+0.1*np.exp(0.05j*2.0*np.pi*x)
#plt.plot(x,y,lw=2,label='sin(x)')
#plt.xlabel('$Time$')
#plt.ylabel('$C(t)$')
#dat0 = ' '.join(map(str, corr.tolist()))
#dat1 = ' '.join(map(str, cori.tolist()))
dat = ''
for i in range(0,len(x)):
    dat = dat+str(x[i])+'+'+str(y[i])+'i ' 
    
#print(dat)
print('time step',t[2]-t[1])
#plt.subplot(2,1,2)
#dat = np.genfromtxt(fname='/home/bing/gwp/spo_2d/1.0.0/cor1') 
f = open('test.dat', 'w')
f.write(str(dat))
#f2 = open('harm_cori.dat','w')
#f2.write(str(dat1))
#np.savetxt('harm_cor.dat',str(dat),fmt="%s")
#for x in range(1,3):
#plt.plot(dat[:,0],dat[:,1],'--',label='$\Re(C(t))$',lw=2)
#plt.plot(dat[:,0],dat[:,2],'--',label='$\Im(C(t))$',lw=2)
#z = np.sqrt(data[:,1]**2+data[:,2]**2)
#plt.plot(data[:,0],z,label='$|C(t)|$',lw=1)
#plt.ylim(-0.2,0.2) 
#plt.legend()
#plt.xlim(0,36)
#plt.savefig('cor.pdf')
#plt.show() 

