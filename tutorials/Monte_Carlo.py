# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from atom import Atom

import matplotlib.pyplot as plt

# we create the system
system = System("2D", K1=0)
system.dipolar = False


# we add some atoms
a = 1.0
N = 10
for i in range(N):
    for j in range(N):
        system.add_object(Atom([i*a, j*a], [0, 1], r=0.2, M=2.2))
        

# we couple every atoms with their neighbors
system.coupleRadius(a+0.1 , 0.1)

# we randomize every magnetizations
system.randomize_magnetizations()

# drawing the initial system
system.draw(name="ferro_init", size=N, center=[(N-1)/2.0, (N-1)/2.0], couples=False)

# setting the system temperature
system.T = 50

# Monte Carlo
for i in range(10):
    system.MonteCarlo(1000)

    system.draw(name="ferro_mc{}".format(i), size=N, center=[(N-1)/2.0, (N-1)/2.0], couples=False)



























