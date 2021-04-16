# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from atom import Atom
from statistics import thermal_at, thermal

import matplotlib.pyplot as plt

# we create the system
system = System("2D", K1=0)
system.dipolar = False


# we add some atoms
a = 1.0
N = 2

for i in range(N):
    for j in range(N):
        system.add_object(Atom([i*a, j*a], [1, 1], r=0.2, M=2.2))

# we couple every atoms with their neighbors
system.coupleRadius(a+0.1 , 0.01)

# we calculate the susceptibility and the specific heat of the system
Ts, xi, cv, m = thermal(system, N=50, Nsample=50, T1=5, T2=300)

plt.plot(Ts, xi)
plt.savefig("xi.png", dpi=200)
plt.clf()

plt.plot(Ts, cv)
plt.savefig("cv.png", dpi=200)
plt.clf()

plt.plot(Ts, m)
plt.savefig("m.png", dpi=200)
plt.clf()






































