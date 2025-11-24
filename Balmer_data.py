import matplotlib.pyplot as plt
import numpy as np

B=[0.87,0.73,1.05,1.12,1.13,1.41,0.96,1.75,1.25,0.86,0.67,0.74,0.65,0.73,1.61,1.28,0.71,0.80,1.30,0.80,1.41,0.90,1.07]

print(len(B))

co='red'
ht='////'
plt.hist(B,range=(0.25,2.5),bins=18,color=co,label='JWST data',histtype='step',hatch=ht)
plt.xlim(0.4,2.4)
plt.xlabel(r"$F_{4200}/F_{3500}$")
plt.legend()

plt.show()
