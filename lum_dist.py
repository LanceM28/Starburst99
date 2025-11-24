import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0=68 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.31)

z1=np.logspace(0,1.5,100)
dist=np.linspace(0,15,10)
d_l0=cosmo.luminosity_distance(z1-1).value/(z1**2)
rad_to_arsec=180.*60*60/np.pi
A_jwst=0.031/rad_to_arsec #JWST
A_elt=A_jwst/12. #extremely large telescope?
fig, (ax2,ax) = plt.subplots(1, 2, figsize=[5.5, 2.8])
#ax.plot(dist,dist*A_jwst*1e6)
#ax.set_xlim([0,10])
#ax.set_ylim([1,200])

ax.plot(z1-1,d_l0*A_elt*1e6/100)
ax2.plot(z1-1,d_l0*A_jwst*1e6/100)

# dist=d_l0/200.0
zz1=11
print('dist(Mpc)=',cosmo.luminosity_distance(zz1-1).value/(zz1**2)/100.0)
#ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.set_xlabel('redshift')
ax.set_xlabel('redshift')
ax2.set_xlim([1,14])
ax2.set_ylim([1,2.7])
ax.set_ylim([0.08,0.23])
ax.set_xlim([1,14])
ax2.set_ylabel('JWST Pixel scale (pc)')
ax.set_ylabel('ELT Pixel scale (pc)')

plt.show()
