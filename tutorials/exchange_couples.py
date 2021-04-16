# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from island import Island

# we create our system
system = System("2D")
    
# we add two islands
system.add_object(Island([-2.5, 0], [-1, 1], a=5.0, b=5.0/1.5, h=2.0, angle=0))
system.add_object(Island([+2.5, 0], [1, -1], a=5.0, b=5.0/1.5, h=2.0, angle=90))
system.add_object(Island([0, 5], [1, 1], a=5.0, b=5.0, h=2.0, angle=0))

# linking the islands
system.couple(0, 2, +0.2)    # +0.2 Ev, ferro
system.couple(0, 1, -0.4)    # -0.4 eV, antiferro

# showing the system caracteristics
print(system.objects[0].caracteristics)
print(system.caracteristics)

# we relax the system and plot it with the links
system.draw(name="couples_init", size=15, center=[0, 2.5], couples=True)

system.relax()

system.draw(name="couples_relaxed", size=15, center=[0, 2.5], couples=True)

