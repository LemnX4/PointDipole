# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import numpy as np

from obj import Object


class Atom(Object):
    def __init__(self, position, magnetization=[1, 0], r=1.0, M=2.2):
        super().__init__(position, magnetization, r, r, r, 0, 0)
        
        self.atomic = True
        self.island = False
        
        self.radius = r
        self.M = M
        self.m = M
        
        self.V = 4/3.0 * np.pi * (r*1e-9)**3
        
        self.Ku = 0.0
        self._u_axis = [0, 0]
        
        
        self.update_values()
        self.update_caracteristics()
        
    @property
    def M(self):
        return self._M
    
    @M.setter
    def M(self, value):
        self._M = value*9.27400949e-24
        self.m = value
        
    @property
    def u_axis(self):
        return self._u_axis
    
    @u_axis.setter
    def u_axis(self, value):
        if np.linalg.norm(value) == 0:
            self.Ku = 0
            self._u_axis = value
        else:
            self._u_axis = value / np.linalg.norm(value)
            





















