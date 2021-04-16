# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8

import numpy as np

from obj import Object


class Island(Object):
    def __init__(self, position, magnetization, a=1.0, b=1.0, h=1.0, angle=0, M0=1720e3):
        super().__init__(position, magnetization, a, b, h, angle, M0)
        
        self.atomic = False
        self.island = True

        self.update_values()
        self.update_caracteristics()