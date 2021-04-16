# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from island import Island

import numpy as np

system = System("2D")

R = 4
N = 6

ar=1.5

system.add_object(Island([0, 0], [1, 0], a=4.0, b=4.0/ar, h=2.0, angle=0))

for i in range(N):
    theta = 2*np.pi/N * i
    
    x = R*np.cos(theta)
    y = R*np.sin(theta)
    
    system.add_object(Island([x, y], [1, 0], a=4.0, b=4.0/ar, h=2.0, angle=0))

for i in range(N):
    theta = 2*np.pi/N * i + np.pi/6
    
    x = np.sqrt(3)*R*np.cos(theta)
    y = np.sqrt(3)*R*np.sin(theta)
    
    system.add_object(Island([x, y], [1, 0], a=4.0, b=4.0/ar, h=2.0, angle=0))
    
for i in range(N):
    theta = 2*np.pi/N * i + np.pi/3
    
    x = 2*R*np.cos(theta)
    y = 2*R*np.sin(theta)
    
    system.add_object(Island([x, y], [1, 0], a=4.0, b=4.0/ar, h=2.0, angle=0))


system.randomize_magnetizations()
system.randomize_angles()

system.couple(0, 1, 5.0)
system.couple(9, 15, -5.0)


system.objects[0].frozen = True


system.draw(name="corona_init", size=20, couples=True, annoted=True)

system.relax()

system.draw(name="corona_relaxed", size=20, couples=True)





























