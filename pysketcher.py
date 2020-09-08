import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

m = 1.  # particle's mass
b = 1
k = 0.0
g = 9.81  # gravity acceleration

# The initial position is (0, 0).
v0 = np.zeros(4)
# The initial speed vector is oriented
# to the top right.
v0[2] = 10.
v0[3] = 0.

def f(v, t0, k):
    # v has four components: v=[u, u'].
    u, udot = v[:2], v[2:]
    # We compute the second derivative u'' of u.
    udotdot = -(b/m)*udot #- (k/m)*u
    #udotdot[1] -= g
    # We return v'=[u', u''].
    return np.r_[udot, udotdot]

fig, ax = plt.subplots(1, 1, figsize=(8, 4))

# We want to evaluate the system on 30 linearly
# spaced times between t=0 and t=3.
t = np.linspace(0., 50., 300)

v = spi.odeint(f, v0, t, args=(k,))
print(v)
ax.plot(v[:, 0], v[:, 1], 'o-', mew=1, ms=6,
            mec='w', label=f'k={k:.1f}')

ax.legend()
ax.set_xlim(-10, 20)
ax.set_ylim(-15, 15)
plt.grid(); plt.show()