# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from island import Island
from kerr import MOKE

# we create the system
system = System("2D")

system.ea1 = [1, 1]
system.ea2 = [1, -1]

R = 7
N = 6

# we add some islands
system.add_object(Island([-3, 0], [1, 0.1], a=5.0, b=4.0, h=2.0, angle=0))
system.add_object(Island([3, 0], [1, 0.1], a=5.0, b=4.0, h=2.0, angle=45.0))
system.add_object(Island([0, 6], [1, 0.1], a=5.0, b=4.0, h=2.0, angle=130.0))


# MOKE measurement of the system and drawing it
MOKE(system, N=100, theta=0.0, Bmax=0.3, plotname="islands_kerr", transverse=True,
     filename="islands_kerr_data", animated=True, anim_size=15, anim_center=[0, 3])

# showing some system caracteristics
print(system.caracteristics)












































