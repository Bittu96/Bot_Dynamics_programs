import math
from math import *
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(8, 4))

m = 5.  # particle's mass
b = 10
k = 5
g = 9.81  # gravity acceleration
c1= 5
c2= 10
lamda1 = (-b + sqrt(b**2 - 4*m*k))/ 2*m
lamda2 = (-b - sqrt(b**2 - 4*m*k))/ 2*m

t = np.linspace(0., 100, 100)
#jz= thrust*jt - g*jt**2 
#jz = (c1*(e**(lamda1*t)) + c2*(e**(lamda2*t)) ) 
lamda1=-0.1
jz = (c1*t*(e**(lamda1*t))) 

ax.plot(t, jz, 'o-', mew=1, ms=6,
            mec='w', label=f'k={k:.1f}')

ax.legend()
ax.set_xlim(-5, 60)
ax.set_ylim(-5, 20)
plt.grid(); plt.show()