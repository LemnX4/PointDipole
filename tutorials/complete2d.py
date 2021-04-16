# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")


from system import System
from island import Island

import numpy as np
import random

# we create the system
system = System("2D")

R = 5
N = 6


# we create a nanowatch
system.add_object(Island([0, 0], [2*random.random()-1, 2*random.random()-1],
                             a=4.0, b=4.0/1.2, h=2.0, angle=random.random()*360))

for i in range(N):
    theta = 2*np.pi/N * i
    
    x = R*np.cos(theta)
    y = R*np.sin(theta)
    
    ar = random.random()*0.5 + 1.0
    
    system.add_object(Island([x, y], [2*random.random()-1, 2*random.random()-1],
                             a=4.0, b=4.0/ar, h=2.0, angle=random.random()*360))

# setting the magnetoctistalline anisotropy
system.K1 = 48e3
system.ea = [1, 2]

# creating some exchange coupling
system.couple(0, 1, 0.02)    # ferro coupling
system.couple(0, 3, -0.01)   # antiferro coupling

# we set a applied field in the 100 direction at -0.1 T
system.Bu = [3, 1]
system.B = -0.1

# we freeze the island number 1
system.objects[1].frozen = True

# we plot our initialized system with couples and annotations
system.draw(name="complete_example_init", size=15, couples=True, annoted=True)

# relaxing it
system.relax()

# we plot the relaxed system
system.draw(name="complete_example_relaxed", size=15, couples=True)

# center island caracteristics
print(system.objects[0].caracteristics)

#system caracteristics
print(system.caracteristics)

































