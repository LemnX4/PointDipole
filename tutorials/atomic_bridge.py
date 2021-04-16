# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from island import Island
from atom import Atom


# we create the system
system = System("3D")
aIsland = 4
aFe = 0.14
N = 10


# we create an atomic bridge with three Fe atoms between two islands
system.add_object(Island([0, 0, 0], [0, 1, 0],
                            a=4.0, b=4.0, h=2.0, angle=0))

for i in range(N):
    pos = aFe*(i+1)
    
    x = pos
    y = 0

    atom = Atom([x+2, y, 0], [1, 0, 0], r=aFe/2.0, M=2.2)
    atom.randomize_magnetization()
    
    system.add_object(atom)

system.add_object(Island([(N+1)*aFe + aIsland, 0, 0], [0, -1, 0],
                            a=4.0, b=4.0, h=2.0, angle=0))


# setting the magnetoctistalline anisotropy to 0
system.K1 = 0

# we do not take into account the dipolar interaction
system.dipolar = False

# creating some exchange coupling
J = 0.01
for i in range(N+1):
    system.couple(i, i+1, J)


# we freeze the nano-islands magnetization
system.objects[0].frozen = True
system.objects[N+1].frozen = True


x = ((N+1)*aFe + aIsland)/2.0
# we plot our initialized system with couples and annotations
system.draw(name="atomic_init", size=x, center=[x, 0], couples=False, annoted=True)

# relaxing it
system.relax(gtol=1e-10)

# we plot the relaxed system
system.draw(name="atomic_relaxed", size=x, center=[x, 0], couples=False)

# atom caracteristics
print(system.objects[1].caracteristics)

#system energies
print(system.energies)

#system caracteristics
print(system.caracteristics)

































