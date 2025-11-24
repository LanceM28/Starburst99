import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def Htime(z):
        t0=566.0
        return t0*((1+z)/10.)**-1.5

def Hredshift(t):
        t0=566.0
        return 10.0*(t/t0)**-(2./3)-1.0

xx=4
# wav = amstrong, y1 = erg/s/Amstrong
# y1= 3  Myr
# y2= 5 Myr
n1=1
n2=2
NN=400
dim=1221
nrows=dim*NN
time,wav, y1 = np.loadtxt('output49291/Kroupa-lowmet1.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1 = np.loadtxt('output15012/Kroupa.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1 = np.loadtxt('output61239/Kroupa-lowmet.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1 = np.loadtxt('output61239/Kroupa-lowmet.spectrum1', unpack=True, usecols=(0,n1,n2))
#time,wav, y1 = np.loadtxt('output38744/top-heavy.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1b = np.loadtxt('output54345/top-heavy-lowmet.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
time=time/1e6
imax=len(y1)/dim
print(imax)

z1,ra,mass1,n_H,met=np.loadtxt('logSFC-CCFid.txt',unpack=True, usecols=(2,4,5,8,9))
z2,ra,mass2,n_H,met=np.loadtxt('logSFC-low1.txt',unpack=True, usecols=(2,4,5,8,9))
z3,ra,mass3,n_H,met=np.loadtxt('logSFC-high1.txt',unpack=True, usecols=(2,4,5,8,9))
Ht1=Htime(z1)
Ht2=Htime(z2)
Ht3=Htime(z3)


fig,axs = plt.subplots(2, 1)
plt.subplots_adjust(hspace=0.6)


ndim=400
for ii in range(1,4):
  if (ii==1):
    z=z1
    Ht=Ht1
    mass=mass1
    co='blue'
    run='CC-Fid'
    ht='////'
  if (ii==2):
    z=z2
    Ht=Ht2
    mass=mass2
    co='green'
    run='SFE 0.35'
    ht='....'
  if (ii==3):
    z=z3
    Ht=Ht3
    mass=mass3
    co='red'
    run='SFE 0.70'
    ht='\\\\'
  htmax=np.max(Ht)-10.0
  htmin=np.min(Ht)+10.0
  z_res=np.zeros(ndim)
  x_res=np.zeros(ndim)
  y_res=np.zeros(ndim)
  for k in range(ndim):
    time0=htmin+np.random.rand()*(htmax-htmin)
    #redsh0=Hredshift(time0)
    redsh0=7.0

    t0=time0-Ht
    tmax0=np.max(t0)
    tmin0=np.min(t0)

    #cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    #d_L0=cosmo.luminosity_distance(redsh0).value*3.08e24

    # in micron
    wav_obs0=wav*(redsh0+1)*1e-4
    #4200 pm 50
    #3500 pm 100
    #3500->i=357
    #4200->i=2055
    #print(wav[356-5],wav[356],wav[356+4])
    #print(wav[392-3],wav[392],wav[392+2])
    #print(wav[182])
    # F=L/(4pi*d_L^2) [erg/s/cm2/Amsr]
    # L_lamb/(1+z) L/(1+z)
    #y1_obs0=np.array(10**(y1-np.log10(4.*np.pi)-2*np.log10(d_L0))/(1+redsh0))
    # convert to micro Jy: 1e-23 cgs

    mag=10.0
    nu=3e10/(wav*1e-8)
    conv=wav/nu
    #y1_obs0=y1_obs0*conv/1e-29*mag*(1+redsh0)**2
    y1_obs0=mag*np.array(10**(y1))*conv
  
    y_tot0=np.zeros(dim-1)
    for ii in range(0,len(t0)-1):
            i=int(t0[ii]/2.0)
            if (i >=0 and i<imax):
                #print(ii,i)
                a1=dim*i
                b1=dim*(i+1)-1
                weight=mass[ii]/1.e6
                y_tot0=y_tot0+weight*y1_obs0[a1:b1]    
                #print(wav_obs[a1],wav_obs[b1])
  
    m1,in1=np.polyfit(wav_obs0[356-5:356+4],y_tot0[356-5:356+4],1)
    B_left=m1*0.35*(redsh0+1)+in1
    m2,in2=np.polyfit(wav_obs0[392-3:392+2],y_tot0[392-3:392+2],1)
    B_right=m2*0.42*(redsh0+1)+in2

    #B_left=np.sum(y_tot0[356-50:356+50])/100.0
    #B_right=np.sum(y_tot0[392-25:392+25])/50.0
    z_res[k]=redsh0
    x_res[k]=y_tot0[182]
    y_res[k]=B_right/B_left
    print('res=',time0, y_tot0[182],y_tot0[392]/y_tot0[356])
    print('res=',B_right/B_left)
  axs[0].scatter(x_res, y_res,color=co,label=run)
  axs[1].hist(y_res,range=(0.5,2.5),bins=16,color=co,label=run,histtype='step',hatch=ht)

axs[0].set_xscale('log')
axs[0].set_yscale('linear')
#axs[0].set_xlim(0.5, 6)
axs[0].set_ylim(0,2.5)
axs[0].set_ylabel(r"$F_{4200}/F_{3500}$")
axs[0].set_xlabel(r"L$_\nu(ergs/s/Hz)$")
axs[0].legend()
axs[1].set_xlabel(r"$F_{4200}/F_{3500}$")
axs[1].legend()

'''
axs[1].plot(wav_obs0[a1:b1], y_tot0, label=f'{run:s} t={time0:3.0f} Myr (z={redsh0:2.2f})')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_ylim(1e-6,1e-3)
axs[1].plot([wav_obs0[392+2],wav_obs0[392+2]],[1e-20,1e-3],c='red')
axs[1].plot([wav_obs0[392-3],wav_obs0[392-3]],[1e-20,1e-3],c='red')
axs[1].plot([wav_obs0[392],wav_obs0[392]],[1e-20,1e-3],c='red')
axs[1].plot([wav_obs0[356-5],wav_obs0[356-5]],[1e-20,1e-3],c='green')
axs[1].plot([wav_obs0[356+4],wav_obs0[356+4]],[1e-20,1e-3],c='green')
axs[1].plot([wav_obs0[356],wav_obs0[356]],[1e-20,1e-3],c='green')
'''

plt.show()
