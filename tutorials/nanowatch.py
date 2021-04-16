# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from island import Island
from kerr import MOKE

import numpy as np

system = System("2D")

R = 5
N = 6

angles = [25, 311, 98, 63, 78, 152, 2, 189]

system.add_object(Island([0, 0], [1, 0], a=4.0, b=4.0, h=2.0, angle=angles[7]))

for i in range(N):
    theta = 2*np.pi/N * i
    
    x = R*np.cos(theta)
    y = R*np.sin(theta)
    
    system.add_object(Island([x, y], [1, 0], a=4.0, b=3.0, h=2.0, angle=angles[i]))


system.draw(name="nanowatch_init", size=15)

system.randomize_magnetizations()

MOKE(system, theta=45.001, Bmax=0.25, N=100, plotname="nanowatch_moke")
    
system.draw(name="nanowatch_relaxed", size=15)

































