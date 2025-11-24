import matplotlib.pyplot as plt
import numpy as np

xx=4
# wav = amstrong, y1 = erg/s/Amstrong
# y1= 3  Myr
# y2= 5 Myr
n1=1
n2=2
dim=1221
nrows=dim*8
time,wav, y1 = np.loadtxt('output49291/Kroupa-lowmet1.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1 = np.loadtxt('output15012/Kroupa.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1 = np.loadtxt('output61239/Kroupa-lowmet.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1 = np.loadtxt('output38744/top-heavy.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
#time,wav, y1b = np.loadtxt('output54345/top-heavy-lowmet.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
time,wav, y1b = np.loadtxt('output49291/Kroupa-lowmet1.spectrum1', unpack=True, usecols=(0,n1,n2),max_rows=nrows)
boost=y1b[dim*xx+dim//5]-y1[dim*xx+dim//5]
y1=y1+boost
time=time/1e6

def Halpha(QHI):
    #erg/s
    # xi_ion=QHI/Lv1500[s-1/(erg s-1 Hz-1] [Hz/erg]
    #typical Lv1500~1e27*(M_*/1e8Msun)
    #typical xi_ion=3e25
    #typical QHI=xi_ion*Lv1500=3e52 s-1(M_*/1e8Msun)
    return 1.36e-12*QHI

def Lv1500(Muv):
    # log10 Lv1500 (erg s−1 Hz−1)
    return 0.4*(51.63-Muv)

def Ll1500(Muv):
    #lam=c/nu
    #L_lam=(c/lam^2)dL/dnu
    # log10 Ll1500 (erg s−1 Ams-1)
    return Lv1500(Muv)+np.log10(3.e10/(1500.*1e-8)/1500.)

Muv=-14
# 1e6 Msun starbust Ll1500 = 39 - 39.5 the first few Myr
print(f"Mag = {Muv:f} -> Solar masses (10Myr): {1e6*10**(Ll1500(Muv)-39):e}")
#EW > 2080 A
#EW ~3000 A for Salpeter IMF to 100 Msun t<5 Myr

mag=120
mass=5e3

xf=mag*mass/1e6
z=6.639
# in cm (Ned's calc)
d_L=66348.7*3.08e24

# in micron
wav_obs=wav*(z+1)*1e-4
# F=L/(4pi*d_L^2) [erg/s/cm2/Amsr]
y1_obs=xf*10**(y1-np.log10(4.*np.pi)-2*np.log10(d_L))

y1b_obs=xf*10**(y1b-np.log10(4.*np.pi)-2*np.log10(d_L))

a1=dim*xx
b1=dim*(xx+1)
plt.plot(wav_obs[a1:b1], y1_obs[a1:b1], label=f"Stellar Kroupa t={time[a1+1]:.1f} Myr")
plt.plot(wav_obs[a1:b1], y1b_obs[a1:b1], label='Stellar Top-heavy')
plt.yscale('log')
plt.xlim(0.5, 6)
plt.ylim(1e-24,1e-18)
plt.xlabel(r"$\lambda_{obs}/\mu$m")
plt.ylabel(r"F$_\lambda(erg/s/cm^2/A)$")
plt.legend()
print(time[dim-1],time[dim*2-1],time[dim*3-1],time[dim*4-1],time[dim*5-1],time[dim*6-1],time[dim*7-1],time[dim*8-1])
plt.show()
