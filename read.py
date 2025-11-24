import matplotlib.pyplot as plt
import numpy as np

def Htime(z):
        t0=566.0
        return t0*((1+z)/10.)**-1.5

#redshift,radius,mass,n_H,met=np.loadtxt('logSFC-low1.txt',unpack=True, usecols=(2,4,5,8,9))
#redshift,radius,mass,n_H,met=np.loadtxt('logSFC-high1.txt',unpack=True, usecols=(2,4,5,8,9))
redshift,radius,mass,n_H,met=np.loadtxt('logSFC-CCFid.txt',unpack=True, usecols=(2,4,5,8,9))
Ht=Htime(redshift)

#SFR
tstart=300
tend=700
nbin=int((tend-tstart)/1.0)
print('nbins=',nbin)
ti=np.linspace(300,700,nbin)
dti=(ti[1:]-ti[:-1])*1e6

massc=np.cumsum(mass)
massci=np.interp(ti,Ht,massc)
sfri=(massci[1:]-massci[:-1])/dti

plt.plot(ti[1:],sfri)
plt.show()
