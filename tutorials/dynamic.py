# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import sys
sys.path.insert(1, "../src")

from system import System
from atom import Atom

import numpy as np
import matplotlib.pyplot as plt


# we create the system
system = System("3D", K1=0)
system.dipolar = False

# we add an atom to the system
system.add_object(Atom([0, 0, 0], [1, 0, 0], r=0.1, M=1.0))

# we set the applied magnetic field and temperature
system.Bu = [0, 0, 1]
system.B = 1

system.T = 0.001    # 1 mK

# LLG
system.LLGevolve(time_span=5e-10, dt=1e-12, gamma=1.7e11, alpha=0.1)

# gathering the time dependent data
magx = []
magxy = []
magz = []
time = []


for i in range(len(system.time)):
    time.append(system.time[i]/1e-12)
    x = system.objects[0].mag_history[i][0]*system.objects[0].M/9.274e-24
    y = system.objects[0].mag_history[i][1]*system.objects[0].M/9.274e-24
    z = system.objects[0].mag_history[i][2]*system.objects[0].M/9.274e-24
    magx.append(x)
    magxy.append(np.hypot(x, y))
    magz.append(z)


# plot the data
plt.plot(time, magx, label="$m_x$")
plt.plot(time, magz, label="$m_z$")

plt.plot(time, magxy, color="C0", linestyle="dashed", alpha=0.25)
plt.axhline(y=0, color="black", linestyle="dotted", alpha=0.25)

plt.xlabel("Time (ps)")
plt.ylabel("Projected atom magnetization ($µ_B$)")

plt.legend()
plt.savefig("time_evolution", dpi=200)












