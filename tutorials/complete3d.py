# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")


from system import System
from island import Island
from kerr import MOKE

import numpy as np
import random

# we create the system
system = System("3D")

R = 5
N = 6


# we create a nanowatch
system.add_object(Island([0, 0, 0], [1, 0, 0], a=4.0, b=4.0/1.2, h=2.0, angle=0.0))

for i in range(N):
    theta = 2*np.pi/N * i
    
    x = R*np.cos(theta)
    y = R*np.sin(theta)
    
    ar = random.random()*0.5 + 1.0
    
    system.add_object(Island([x, y, 0], [1, 0, 0], a=4.0, b=4.0/ar, h=2.0, angle=0.0))

system.randomize_magnetizations()
system.randomize_angles()


# setting the magnetoctistalline anisotropy
system.K1 = 48e3

# creating some exchange coupling
system.couple(0, 1, 0.002)    # ferro coupling
system.couple(0, 3, -0.001)   # antiferro coupling

# we set a applied field in the 311 direction at -0.1 T
system.Bu = [3, 1, 1]
system.B = -0.1

# we freeze the island number 1
system.objects[1].frozen = True


# we plot our initialized system with couples and annotations
system.draw(name="complete3d_init", size=15, couples=True, annoted=True)

# relaxing it
system.relax()

# we plot the relaxed system
system.draw(name="complete3d_relaxed", size=15, couples=True)


# center island caracteristics
print(system.objects[0].caracteristics)

#system caracteristics
print(system.caracteristics)

































