# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import numpy as np
import random

from demag import demagnetization

from llg import Bth

class Object:
    def __init__(self, position, magnetization, a=1.0, b=1.0, h=1.0, angle=0, M0=1720e3):
        if len(magnetization) == 2 and len(position) == 2:
            self.dim = "2D"
        elif len(magnetization) == 3 and len(position) == 3:
            self.dim = "3D"
        elif len(magnetization) != len(position):
            print("\nError : island dimension undefined.\n")
            self.dim == "undefined"
        
        self._frozen = False
        self._atomic = False
        self._island = False
        
        self.pos = position
        
        self._a = a*1e-9
        self._b = b*1e-9
        self._h = h*1e-9
        
        self.mag = magnetization
        
        self.thermic_field = []
        self.time_evolved = False
        self.time = []
        self.mag_history = []
        
        self.angle = angle
        self.M0 = M0
        
        self.Ku = 0.0
        
        self.har = a/b
        self.var = a/h
        
        self.coupled_with = []
    
    @property
    def frozen(self):
        return self._frozen
    
    @frozen.setter
    def frozen(self, value):
        self._frozen = value
    
    @property
    def atomic(self):
        return self._atomic
    
    @atomic.setter
    def atomic(self, value):
        self._atomic = value
    
    @property
    def island(self):
        return self._island
    
    @island.setter
    def island(self, value):
        self._island = value

    @property
    def a(self):
        return self._a
    
    @a.setter
    def a(self, value):
        self._a = value*1e-9
        self.update_values()
    
    @property
    def b(self):
        return self._b
    
    @b.setter
    def b(self, value):
        self._b = value*1e-9
        self.update_values()
        
    @property
    def h(self):
        return self._h
    
    @h.setter
    def h(self, value):
        self._h = value*1e-9
        self.update_values()
    
    @property
    def M0(self):
        return self._M0
    
    @M0.setter
    def M0(self, value):
        self._M0 = value
        self.update_values()
    
    @h.setter
    def h(self, value):
        self._h = value*1e-9    
        self.update_values()
    
    @property
    def mag(self):
        return self._mag
    
    @mag.setter
    def mag(self, value):
        if np.linalg.norm(value) == 0:
            if self.dim == "2D":
                print("Error : null vector magnetization. Set to [1, 0] by default.")
                self._mag = [1, 0]
            elif self.dim == "3D":
                print("Error : null vector magnetization. Set to [1, 0, 0] by default.")
                self._mag = [1, 0, 0]
        else:
            self._mag = value / np.linalg.norm(value)
    
    def randomize_magnetization(self, angle=180):
        theta = 2*(0.5-random.random())*angle*np.pi/180
        if self.dim == "2D":
            new_angle = np.arctan2(self.mag[1], self.mag[0]) + theta
            self.mag = [np.cos(new_angle), np.sin(new_angle)]
        elif self.dim == "3D":
            self.mag = [random.random()-0.5, random.random()-0.5, random.random()-0.5]
    
    def update_values(self):
        if self.atomic:
            self.factors = [1/3.0, 1/3.0, 1/3.0]
            self.ku = 0
            return
        
        self.V = np.pi * (self.a/2.0)*(self.b/2.0)*self.h
        self.factors = demagnetization(self.a, self.b, self.h)
        self.M = self.V*self.M0
        
        if self.dim == "2D":
            self.E0 = (2*np.pi*1e-7) * self.M0**2 * self.V * self.factors[0]
            self.ku = self.M0**2 *(2*np.pi*1e-7)*(1 - 2*self.factors[0] - self.factors[2])

    
    def update_caracteristics(self):
        if self.atomic:
            c = "\n##########\tCaracteristics of the atom:\t##########\n\n"
        elif self.island:
            c = "\n##########\tCaracteristics of the nano-island:\t##########\n\n"
            
        c += "Position (nm) : {}\n".format(self.pos)
        if self.atomic:
            c += "Magnetic moment : {} µB\n".format(self.m)
        else:
            c += "Magnetic moment : {} µB\n".format(self.M/9.74e-24)
        
        f = ""
        if self.frozen:
            f = " (frozen)"
        c += "Magnetization direction{} : {}\n".format(f, self.mag)
        
        if not self.atomic:
            c += "Angle : {}°\n".format(self.angle)
            c += "Large diamater (a) : {} nm\n".format(self.a/1e-9)
            c += "Small diamater (b) : {} nm\n".format(self.b/1e-9)
            c += "Height (h) : {} nm\n".format(self.h/1e-9)
            c += "Horizontal aspect ratio : {} \n".format(self.har)
            c += "Vertical aspect ratio : {} \n".format(self.var)
        else:
            c += "Radius : {} nm\n".format(self.radius)
            
        c += "Volume : {} m³\n".format(self.V)
        
        if not self.atomic:
            if self.dim == "2D":
                c += "Uniaxial horizontal constant Ku : {} kJ/m³\n".format(self.Ku/1e3)
                c += "Self dipolar energy : {} J\n".format(self.E0)
            elif self.dim == "3D":
                c += "Demagnetization factors : {}\n".format(self.factors)
        elif self.Ku != 0:
            c += "Uniaxial constant Ku : {} kJ/m³\n".format(self.Ku/1e3)
            c += "Uniaxial axis : {}\n".format(self.u_axis)
            
        if len(self.coupled_with) !=0 :
            c += "\nCoupled with islands : {}\n".format(self.coupled_with)
        else:
            c += "\nNot coupled with other islands.\n"
            
        c += "\n###################################################################\n"
        
        self.caracteristics = c

    def initialize_thermic_field(self, gamma, alpha, T, N):
        self.thermic_field = []
        for i in range(N+1):
            self.thermic_field.append(Bth(gamma, alpha, self.M, T))




















