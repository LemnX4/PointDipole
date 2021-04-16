# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import random
from scipy import optimize
from scipy.integrate import solve_ivp
import copy

from exchange import Couple
from energy import *
from dipolar import *

from llg import LLG

class System:
    def __init__(self, dimension, K1=48e3):
        self.dim = dimension
        
        if self.dim == "2D":
            self._ea1 = [1, 0]
            self._ea2 = [0, 1]
            
            self._Bu = [0, 0]
            
        elif self.dim == "3D":
            self._ea1 = [1, 0, 0]
            self._ea2 = [0, 1, 0]
            self._ea3 = [0, 0, 1]
            
            self._Bu = [0, 0, 0]
            
        else:
            print("\nError: the system should have a 2D or 3D dimension.\n")
        
        self.E_biaxial = 0.0
        self.E_uniaxial = 0.0
        self.E_shape = 0.0
        self.E_dipolar_interaction = 0.0
        self.E_dipolar_self = 0.0
        self.E_dipolar = 0.0
        self.E_Zeeman = 0.0
        self.E_exchange = 0.0
        self.E_total = 0.0
        
        self.B = 0.0
        self.K1 = K1
        
        self.T = 0.0
        self.V = 0.0
        
        self.time_evolved = False
        self.time = []
        
        self.objects = []
        self.couples = []
        self.bridges = []
        
        self.dipolar = True

    @property
    def ea1(self):
        return self._ea1

    @ea1.setter
    def ea1(self, value):
        if np.linalg.norm(value) == 0:
            if self.dim == "2D":
                print("Error : null vector easy axis. Set to [1, 0] by default.")
                self._ea1 = [1, 0]
            else:
                print("Error : null vector easy axis. Set to [1, 0, 0] by default.")
                self._ea1 = [1, 0, 0]
        else:
            self._ea1 = value/np.linalg.norm(value)
            
        self.update_caracteristics()
        
    @property
    def ea2(self):
        return self._ea2

    @ea2.setter
    def ea2(self, value):
        if np.linalg.norm(value) == 0:
            if self.dim == "2D":
                print("Error : null vector easy axis. Set to [1, 0] by default.")
                self._ea2 = [1, 0]
            else:
                print("Error : null vector easy axis. Set to [1, 0, 0] by default.")
                self._ea2 = [1, 0, 0]
        else:
            self._ea2 = value/np.linalg.norm(value)
            
        self.update_caracteristics()
            
    @property
    def ea3(self):
        return self._ea3

    @ea3.setter
    def ea3(self, value):
        if np.linalg.norm(value) == 0:
            if self.dim == "2D":
                print("Error : null vector easy axis. Set to [1, 0] by default.")
                self._ea3 = [1, 0]
            else:
                print("Error : null vector easy axis. Set to [1, 0, 0] by default.")
                self._ea3 = [1, 0, 0]
        else:
            self._ea3 = value/np.linalg.norm(value)
            
        self.update_caracteristics()
        
    @property
    def Bu(self):
        return self._Bu
    
    @Bu.setter
    def Bu(self, value):
        if np.linalg.norm(value) != 0:
            self._Bu = value/np.linalg.norm(value)
        else:
            self._Bu = value
        self.update_caracteristics()


    def add_object(self, obj):
        self.objects.append(obj)
        self.update_caracteristics()
        self.update_energies()
        
    def remove_object(self, index):
        self.objects.remove(index)
        self.update_caracteristics()
        self.update_energies()
    
    
    def clear(self, show=True):
        self.objects = []
        self.couples = []
        
        self.update_caracteristics()
        self.update_energies()
        
        if show:
            print("\nSystem cleared.\n")
    
    def randomize_positions(self, delta=0.01):
        for obj in self.objects:
            for i in range(len(obj.pos)):
                obj.pos[i] += (0.5-random.random())*delta 
                
    def randomize_magnetizations(self, angle=180):
        for obj in self.objects:
            obj.randomize_magnetization(angle)
        self.update_caracteristics()
            
    def randomize_angles(self):
        for obj in self.objects:
            obj.angle = random.random()*360
        self.update_caracteristics()
    
    def couple(self, index1, index2, J=0.05):
        if J==0:
            return
        for i in range(len(self.couples)):
            i1 = self.couples[i].index1 
            i2 = self.couples[i].index2
            t1 = i1 == index1 and i2 == index2 
            t2 = i1 == index2 and i2 == index1
            if t1 or t2:
                if self.couples[i].J != J:
                    print("\nCouple already exists : J changed from {} eV to {} eV.\n".format(self.couples[i].J, J))
                    self.couples[i].J = J
                #else:
                    #print("\nCouple already exists with the same J.\n")
                return False
            
        self.objects[index1].coupled_with.append(index2)
        self.objects[index2].coupled_with.append(index1)
        
        self.objects[index1].update_caracteristics()
        self.objects[index2].update_caracteristics()
        
        self.couples.append(Couple(index1, index2, J))
        self.update_caracteristics()
        return True
    
    def uncouple(self, object_index):
        self.objects[object_index].coupled_with = []
        self.objects[object_index].update_caracteristics()
        
        to_remove = []
        for i in range(len(self.bridges)):
            couple = self.couples[i]
            
            remove = False
            if couple.index1 == object_index:
                self.objects[couple.index2].coupled_with.remove(couple.index1)
                self.objects[couple.index2].update_caracteristics()
                remove = True
                
            if couple.index2 == object_index:
                self.objects[couple.index1].coupled_with.remove(couple.index2)
                self.objects[couple.index1].update_caracteristics()
                remove = True
            
            if remove:
                to_remove.append(couple)
        
        new_couples = []
        for b in self.couples:
            if b not in to_remove:
                new_couples.append(b)
        
        self.couples = new_couples
        print("\nObject {} uncoupled.\n".format(object_index))
        self.update_caracteristics()
    
    
    def remove_couple(self, index1, index2):
        for i in range(len(self.couples)):
            couple = self.couples[i]
            
            if couple.index1 == index1 and couple.index2 == index2:
                self.objects[couple.index1].coupled_with.remove(couple.index2)
                self.objects[couple.index2].coupled_with.remove(couple.index1)
            
            self.couples.remove(couple)
            print("\nCouple between object {} and {} removed.\n".format(index1, index2))
            self.update_caracteristics()
            return
        print("\nNo couple to remove.\n")
    
    def coupleRadius(self, radius, J=0.1):
        for i in range(len(self.objects)):
            for j in range(len(self.objects)):
                if i==j:
                    continue
                r = np.array(self.objects[i].pos) - np.array(self.objects[j].pos)
                d = np.linalg.norm(r)
                
                if d < radius:
                    self.couple(i, j, J)
                    
        
    def update_caracteristics(self):
        c = "\n##########\tCaracteristics of the system:\t##########\n\n"
        c += "Number of nano-objects : {}\n".format(len(self.objects))
        c += "Magnetocristalline anisotropy constant : {} kJ/m³\n".format(self.K1/1e3)
        
        c += "Easy axis 1 : {}\n".format(self.ea1)
        
        if self.dim == "2D":
            c += "Easy axis 2 : {}\n\n".format(self.ea2)
        elif self.dim == "3D":
            c += "Easy axis 2 : {}\n".format(self.ea2)
            c += "Easy axis 3 : {}\n\n".format(self.ea3)
            
        if self.B == 0 or np.linalg.norm(self.Bu) == 0:
            c += "No applied field.\n\n".format(self.B)
        else:
            c += "Applied field axis : {}\n".format(self.Bu)
            c += "Applied field value : {} T\n\n".format(self.B)
        
        if self.T != 0:
            c += "System temperature : {} K\n\n".format(self.T)
        
        if len(self.couples) != 0:
            for i in range(len(self.couples)):
                string = "Couple {} : object {} and {} coupled with J = {} eV\n"
                c += string.format(i, self.couples[i].index1, self.couples[i].index2, self.couples[i].J)
            c += "\n"
        else:
            c += "There is no exchange couples.\n\n"
        
        if self.dim == "2D":
            mag = np.array([0.0, 0.0])
        elif self.dim == "3D":
            mag = np.array([0.0, 0.0, 0.0])
            
        self.V = 0.0
        for i in range(len(self.objects)):
            self.V += self.objects[i].V
        
            
        mag_value = 0
        
        for i in range(len(self.objects)):
            mag += np.array(self.objects[i].mag) * self.objects[i].M
        
        self.mag = mag
        
        
        mag_value = np.linalg.norm(mag)
        if  mag_value != 0:
            mag /= mag_value
        
        self.M = mag_value
        
        
        mag_sat = 0
        for i in range(len(self.objects)):
            mag_sat += self.objects[i].M
        
        self.Ms = mag_sat
        
        
        c += "Total magnetization vector (unitary) : {}\n".format(mag)
        c += "Magnetization : {} Am²\n".format(mag_value)
        c += "Saturation magnetization : {} Am²\n\n".format(mag_sat)
        
        c += "Total energy : {} J\n".format(self.E_total)
        
        c += "\n##########################################################\n"
        
        self.caracteristics = c
        return c
        
        
    def draw(self, name="system", scale=1.0, size=4.0, center=[0, 0], fontscale=1.0, display=True,
             anisotropy=True, applied=True, annoted=False, couples=True):
        
        fontsize = 7.0*fontscale
        
        fig = plt.figure(0)
        ax = fig.add_subplot(111, aspect='equal')
        
        
        plt.grid(alpha=0.25)
        
        plt.xlim(center[0]-size/2.0, center[0]+size/2.0)
        plt.ylim(center[1]-size/2.0, center[1]+size/2.0)
        
        scale = scale*size/20.0
        
        # objects
        for i in range(len(self.objects)):
            obj = self.objects[i]
            
            x = obj.pos[0] - obj.mag[0]*scale
            y = obj.pos[1] - obj.mag[1]*scale
            
            plt.arrow(x, y, obj.mag[0]*2*scale, obj.mag[1]*2*scale,
                      width=0.1*scale, head_width=0.3*scale, head_length=0.3*scale,
                      fc='black', ec='black', length_includes_head = True)
            
            front_color = "C0"
            edge_color = "blue"
            
            if couples and len(obj.coupled_with) != 0:
                front_color = "purple"
                edge_color = "darkviolet"
                
            
            if obj.frozen:
                front_color = "slategrey"
                edge_color = "black"
            
            ellipse = Ellipse(obj.pos, width=obj.a/1e-9, height=obj.b/1e-9, angle=obj.angle,
                              fill=False, ec=edge_color)
            ax.add_artist(ellipse)
            
            ellipse = Ellipse(obj.pos, width=obj.a/1e-9, height=obj.b/1e-9, angle=obj.angle,
                              fill=True, alpha=0.2, fc=front_color)
            ax.add_artist(ellipse)
            
            if annoted:
                displacement = 0.6*size/20.0
                projmag = [obj.mag[0], obj.mag[1]]
                projmag /= np.linalg.norm(projmag)
                xt = obj.pos[0] - displacement*projmag[1]
                yt = obj.pos[1] + displacement*projmag[0]
                plt.text(xt, yt, str(i), ha="center", va="center", fontsize=fontsize)
                
                if obj.frozen:
                    xt = obj.pos[0] + displacement*obj.mag[1]
                    yt = obj.pos[1] - displacement*obj.mag[0]
                    plt.text(xt, yt, "F", ha="center", va="center", fontsize=fontsize)
        
        ascale = size/10.0
        
        # exchange couples
        if couples:
            for i in range(len(self.couples)):
                if self.couples[i].J == 0:
                    continue
                
                p1 = self.objects[self.couples[i].index1].pos
                p2 = self.objects[self.couples[i].index2].pos
                
                x = p1[0]
                y = p1[1]
                vx = p2[0]-p1[0]
                vy = p2[1]-p1[1]
                
                theta = np.arctan2(vy, vx)
                a = size/15.0
                
                plt.arrow(x+a*np.cos(theta), y+a*np.sin(theta), vx-2*a*np.cos(theta), vy-2*a*np.sin(theta),
                          width=0.25*ascale, head_width=0, head_length=0,
                          fc='purple', ec='purple', alpha=0.2, length_includes_head = True, shape="full")
                
                l = np.hypot(abs(vx), abs(vy))
                xt = x + l/2.0*np.cos(theta)
                yt = y + l/2.0*np.sin(theta)
                
                if self.couples[i].J > 0:
                    plt.text(xt, yt, "+", ha="center", va="center", color="purple", fontsize=fontsize*1.5)
                else:
                    plt.text(xt, yt, "-", ha="center", va="center", color="purple", fontsize=fontsize*2.5)
        
        # magnetocristalline anisotropy
        if anisotropy and self.K1 != 0:
            x = center[0]+size/2.0-ascale
            y = center[1]-size/2.0+ascale
            
            alpha = 1.0
            
            color = "darkred"
            
            plt.arrow(x, y, 
                         self.ea1[0]*0.8*ascale, self.ea1[1]*0.8*ascale,
                         width=0.03*ascale, head_width=0.2*ascale, head_length=0.2*ascale, alpha=alpha,
                         fc=color, ec=color, length_includes_head = True, shape="full")
            
            plt.arrow(x, y, 
                         -self.ea1[0]*0.8*ascale, -self.ea1[1]*0.8*ascale,
                         width=0.03*ascale, head_width=0.2*ascale, head_length=0.2*ascale, alpha=alpha,
                         fc=color, ec=color, length_includes_head = True, shape="full")
            
            plt.arrow(x, y,
                         self.ea2[0]*0.8*ascale, self.ea2[1]*0.8*ascale,
                         width=0.03*ascale, head_width=0.2*ascale, head_length=0.2*ascale, alpha=alpha,
                         fc=color, ec=color, length_includes_head = True, shape="full")
            
            plt.arrow(x, y,
                         -self.ea2[0]*0.8*ascale, -self.ea2[1]*0.8*ascale,
                         width=0.03*ascale, head_width=0.2*ascale, head_length=0.2*ascale, alpha=alpha,
                         fc=color, ec=color, length_includes_head = True, shape="full")
            
            circle = Circle((x, y), radius=0.1*ascale, ec=color, fc="grey")
            ax.add_artist(circle)
            
        # applied field
        if applied and self.B != 0:
            x = center[0]+size/2.0-ascale
            y = center[1]+size/2.0-ascale
            
            vx = self.Bu[0]*0.8*ascale*np.sign(self.B)
            vy = self.Bu[1]*0.8*ascale*np.sign(self.B)
            
            alpha = 1.0
            
            color = "darkblue"
            
            plt.arrow(x, y, vx, vy,
                         width=0.03*ascale, head_width=0.2*ascale, head_length=0.2*ascale, alpha=alpha,
                         fc=color, ec=color, length_includes_head = True, shape="full")
            
            circle = Circle((x, y), radius=0.1*ascale, ec=color, fc="grey")
            ax.add_artist(circle)
        
        # save
        name += ".png"
        plt.savefig(name, dpi=200)
        if display:
            print('\nSystem image saved as "' + name + '".\n')
        plt.clf()
        return name
    
    def get_magnetizations(self):
        mags = []
        for i in range(len(self.objects)):
            mags.append(self.objects[i].mag)
        return mags
        
    def set_magnetizations(self, mags):
        for i in range(len(self.objects)):
            self.objects[i].mag = mags[i]
    
    def get_magnetization_angles2D(self):
        angles = []
        
        for i in range(len(self.objects)):
            angles.append(np.arctan2(self.objects[i].mag[1], self.objects[i].mag[0]))
            
        return angles
    
    def set_magnetization_angles2D(self, angles):
        for i in range(len(self.objects)):
            self.objects[i].mag = [np.cos(angles[i]), np.sin(angles[i])]
            self.objects[i].update_caracteristics()

    def get_magnetization_angles3D(self):
        angles = []
        
        for i in range(len(self.objects)):
            theta = np.arccos(self.objects[i].mag[2])
            phi = np.arctan2(self.objects[i].mag[1], self.objects[i].mag[0])
            angles.append(theta)
            angles.append(phi)
            
        return angles
    
    def set_magnetization_angles3D(self, angles):
        for i in range(len(self.objects)):
            theta = angles[2*i]
            phi = angles[2*i+1]
            
            self.objects[i].mag[0] = np.sin(theta)*np.cos(phi)
            self.objects[i].mag[1] = np.sin(theta)*np.sin(phi)
            self.objects[i].mag[2] = np.cos(theta)
            
            self.objects[i].update_caracteristics()

    
    def update_energies(self):
        if self.dim == "2D":
            angles = self.get_magnetization_angles2D()
        
            self.E_biaxial = 0.0
            self.E_uniaxial = 0.0
            
            self.E_shape = 0.0
            self.E_dipolar_interaction = 0.0
            self.E_dipolar_self = 0.0
            self.E_dipolar = 0.0
            
            self.E_Zeeman = 0.0
            
            self.E_exchange = 0.0
            
            self.E_total = 0.0
            
            
            for i in range(len(self.objects)):
                mag = [np.cos(angles[i]), np.sin(angles[i])]
                
                # uniaxial
                if self.objects[i].atomic and self.objects[i].Ku != 0:
                    Ku = self.objects[i].Ku
                    if Ku > 0: 
                        self.E_uniaxial += Ku * self.objects[i].V * (1-(np.dot(mag, self.objects[i].u_axis))**2)
                    else:
                        self.E_uniaxial += -Ku * self.objects[i].V *(np.dot(mag, self.objects[i].u_axis))**2
                
                # biaxial
                dots = np.dot(mag, self.ea1)**2 * np.dot(mag, self.ea2)**2
                self.E_biaxial += self.objects[i].V * self.K1 * dots
                
                
                # shape
                if not self.objects[i].atomic and self.dipolar:
                    theta = angles[i] - self.objects[i].angle * np.pi/180
                    self.E_shape += self.objects[i].V * self.objects[i].ku * np.sin(theta)**2
                
                # dipolar
                if self.dipolar:
                    for j in range(i):
                        self.E_dipolar_interaction += dipolar_int2D(self.objects[i], self.objects[j], 
                                                                    angles[i], angles[j]) / 2.0
                    
                    self.E_dipolar_self += self.objects[i].E0
            
                # Zeeman
                m = self.objects[i].M*np.array([np.cos(angles[i]), np.sin(angles[i])])
                B = self.B*np.array(self.Bu)
                self.E_Zeeman += -np.dot(m, B)
        
        
            # exchange
            for i in range(len(self.couples)):
                angle1 = angles[self.couples[i].index1]
                angle2 = angles[self.couples[i].index2]
                mag1 = [np.cos(angle1), np.sin(angle1)]
                mag2 = [np.cos(angle2), np.sin(angle2)]
                self.E_exchange += -self.couples[i].J*np.dot(mag1, mag2)*1.602e-19 
            
            self.E_dipolar = self.E_dipolar_self + self.E_shape + self.E_dipolar_interaction
            
            self.E_total = self.E_dipolar + self.E_biaxial + self.E_uniaxial + self.E_Zeeman + self.E_exchange
            
            e = "\n\tEnergy of the system:\n\n"
            
            e += "Dipolar internal energy : {} J\n".format(self.E_dipolar_self)
            e += "Dipolar uniaxial shape energy : {} J\n".format(self.E_shape)
            e += "Dipolar interaction energy : {} J\n".format(self.E_dipolar_interaction)
            e += "Total dipolar energy : {} J\n\n".format(self.E_dipolar)
            
            if self.E_uniaxial != 0:
                e += "Uniaxial energy : {} J\n\n".format(self.E_uniaxial)
            
            e += "Magnetocristalline energy : {} J\n\n".format(self.E_biaxial)
            
            e += "Zeeman energy : {} J\n\n".format(self.E_Zeeman)
            
            e += "Heisenberg exchange energy : {} J\n\n".format(self.E_exchange)
            
            e += "Total : {} J\n".format(self.E_total)
            
            self.energies = e
            return e
        else:
            angles = self.get_magnetization_angles3D()
        
            self.E_cubic = 0.0
            self.E_uniaxial = 0.0
            
            self.E_shape = 0.0
            self.E_dipolar_interaction = 0.0
            self.E_dipolar_self = 0.0
            self.E_dipolar = 0.0
            
            self.E_Zeeman = 0.0
            
            self.E_exchange = 0.0
            
            self.E_total = 0.0
            
            
            for i in range(len(self.objects)):
                theta = angles[2*i]
                phi = angles[2*i+1]
                
                mag = [np.sin(theta)*np.cos(phi),
                       np.sin(theta)*np.sin(phi),
                       np.cos(theta)]
                
                # uniaxial
                if self.objects[i].atomic and self.objects[i].Ku != 0:
                    Ku = self.objects[i].Ku
                    if Ku > 0: 
                        self.E_uniaxial += Ku * self.objects[i].V * (1-(np.dot(mag, self.objects[i].u_axis))**2)
                    else:
                        self.E_uniaxial += -Ku * self.objects[i].V *(np.dot(mag, self.objects[i].u_axis))**2
                        
                # cubic
                a = np.dot(mag, self.ea1)
                b = np.dot(mag, self.ea2)
                c = np.dot(mag, self.ea3)
                dots = (a*b)**2 + (b*c)**2 + (c*a)**2
                self.E_cubic += self.objects[i].V * self.K1 * dots
                
                
                # shape
                if not self.objects[i].atomic and self.dipolar:
                    phi = self.objects[i].angle*np.pi/180
                    mag_eff = [mag[0]+np.cos(phi), mag[1]-np.sin(phi), mag[2]]
                    mag_eff /= np.linalg.norm(mag_eff)
                    
                    demag = mag_eff[0]**2*self.objects[i].factors[0]
                    demag += mag_eff[1]**2*self.objects[i].factors[1]
                    demag +=  mag_eff[2]**2*self.objects[i].factors[2]
                    self.E_shape += self.objects[i].M0**2 * self.objects[i].V *(2*np.pi*1e-7)*demag
                
                # dipolar
                if self.dipolar:
                    for j in range(i):
                        
                        theta1 = angles[2*i]
                        phi1 = angles[2*i+1]
                        
                        theta2 = angles[2*j]
                        phi2 = angles[2*j+1]
                        
                        self.E_dipolar_interaction += dipolar_int3D(self.objects[i], self.objects[j], 
                                                                    theta1, phi1, theta2, phi2) / 2.0
            
                # Zeeman
                m = self.objects[i].M * self.objects[i].mag
                B = self.B*np.array(self.Bu)
                self.E_Zeeman += -np.dot(m, B)
        
        
            # exchange
            for i in range(len(self.couples)):
                mag1 = self.objects[self.couples[i].index1].mag
                mag2 = self.objects[self.couples[i].index2].mag
                self.E_exchange += -self.couples[i].J*np.dot(mag1, mag2)*1.602e-19 
            
            self.E_dipolar = self.E_shape + self.E_dipolar_interaction
            
            self.E_total = self.E_dipolar + self.E_cubic + self.E_uniaxial + self.E_Zeeman + self.E_exchange 
            
            e = "\n\tEnergy of the system:\n\n"
            
            if self.dipolar:
                e += "Dipolar uniaxial shape energy : {} J\n".format(self.E_shape)
                e += "Dipolar interaction energy : {} J\n".format(self.E_dipolar_interaction)
                e += "Total dipolar energy : {} J\n\n".format(self.E_dipolar)
            else:
                e += "No dipolar energy taken into account.\n\n".format(self.E_dipolar)
            
            if self.E_uniaxial != 0:
                e += "Uniaxial energy : {} J\n\n".format(self.E_uniaxial)
            
            e += "Magnetocristalline energy : {} J\n\n".format(self.E_cubic)
            
            e += "Zeeman energy : {} J\n\n".format(self.E_Zeeman)
            
            e += "Heisenberg exchange energy : {} J\n\n".format(self.E_exchange)
            
            e += "Total : {} J\n".format(self.E_total)
            
            self.energies = e
            return e
    
    
    def relax(self, display=True, gtol=1e-10, epsilon=1e-10):
        
        if display:
            print("\nMinimizing the energy using CG...\n")
        
        if self.dim == "2D":
            angles = self.get_magnetization_angles2D()
        
            new_angles = optimize.fmin_cg(total_energy_min2D, angles, args=(self,), disp=display,
                                          gtol=gtol, epsilon=epsilon)
            
            self.set_magnetization_angles2D(new_angles)
        if self.dim == "3D":
            angles = self.get_magnetization_angles3D()
        
            new_angles = optimize.fmin_cg(total_energy_min3D, angles, args=(self,), disp=display,
                                          gtol=gtol, epsilon=epsilon)
            
            self.set_magnetization_angles3D(new_angles)
            
        if display:
            print("\nEnergy minimized.")
        
        self.update_energies()
        self.update_caracteristics()
        
        if display:
            print(self.energies)
        return

    
    def find_min(self, N=5, detailed=True, display=False):
        N = int(abs(N))
        
        energies = []
        mags = []
        
        if detailed:
            print("Finding minimum for {} different random initial magnetizations...".format(N))

        for i in range(N):
            self.randomize_magnetizations()
            self.relax(display=display)
            energies.append(self.E_total)
            mags.append(self.get_magnetizations())
        
        E_min = min(energies)
        index = energies.index(E_min)
        
        mag_min = mags[index]
        
        self.set_magnetizations(mag_min)
        self.update_caracteristics()
        self.update_energies()
        
        if detailed:
            print("Minimum finding process finished.")
            #print("For {} minimizations, those energies were found (J) :".format(N))
            #print("{}".format(energies))
            print("The minimum energy found is {} J\n".format(E_min))
        
        return E_min
    
    def measure_remanent(self, theta=0.0, Hmax=1.0, Hmin=0.0, N=10, detailed=False, display=True):
        Hdir = [np.cos(theta*np.pi/180), np.sin(theta*np.pi/180)]
        self.Hu = Hdir
        
        N = int(abs(N))
        if N < 2:
            N = 2
        
        Hs = (N-np.linspace(0, N, N))/N*(Hmax-Hmin) + Hmin
        
        for H in Hs:
            self.H = H
            self.relax(display=detailed)
        
        reml = np.dot(self.mag, Hdir)*self.M/self.Ms
        
        Htdir = [Hdir[1], -Hdir[0]]
        remt = np.dot(self.mag, Htdir)*self.M/self.Ms
        
        if display:
            print("\nRemanent magnetization measured for \u03B8={}° :".format(theta, reml, remt))
            print("Longitidinal : {} M$_s$\nTransverse : {} M$_s$\n".format(reml, remt))
        
        return reml, remt

    def LLGevolve(self, time_span=1e-9, dt=1e-12, gamma=2.2e-5, alpha=0.05, rtol=1e-7, display=True):
        if self.dim != "3D":
            print("\nError: 3D is required for LLG time evolution.\n")
            return
        
        if display:
            print("\nApplying the Landau-Lifshitz-Gilbert equation on the system...")
        
        mags = []
        for i in range(len(self.objects)):
            mags.append(self.objects[i].mag[0])
            mags.append(self.objects[i].mag[1])
            mags.append(self.objects[i].mag[2])
        
        if self.T != 0:
            for i in range(len(self.objects)):
                self.objects[i].initialize_thermic_field(gamma, alpha, self.T, int(time_span/dt))
                
        
        result = solve_ivp(LLG, (0, time_span), copy.deepcopy(mags), args=(self, gamma, alpha, dt, 1e-7), method="RK45",
                           first_step=dt, max_step=dt, rtol=rtol)
        
        t = result.t
        mags_t = result.y
        n = len(result.t)-1 
        
        N = int(len(t)/(time_span/dt))
        
        mags = []
        for i in range(len(self.objects)):
            self.objects[i].mag = [mags_t[i*3][n], mags_t[i*3+1][n], mags_t[i*3+2][n]]
            
            for k in range(len(t)):
                if k%N == 0:
                    mag = [mags_t[i*3][k], mags_t[i*3+1][k], mags_t[i*3+2][k]]
                    self.objects[i].mag_history.append(mag)
                    self.objects[i].time.append(t[k])
                    
            self.objects[i].time_evolved = True
                
        for k in range(len(t)):
            if k%N == 0:
                self.time.append(t[k])
                
        self.time_evolved = True
        
        if display:
            print("Finished.\n")

        self.update_energies()
        self.update_caracteristics()

    def MonteCarlo(self, N, display=True):
        old_energy = self.E_total
        kB = 1.380649e-23
        
        if display:
            print("\nDoing Monte Carlo at T={}K...".format(self.T))
        
        for i in range(N):
            index = np.random.randint(len(self.objects))
            if self.objects[index].frozen:
                continue
            old_mag = copy.copy(self.objects[index].mag)
            self.objects[index].randomize_magnetization(angle=180)

            self.update_energies()
            self.update_caracteristics()
            
            current_energy = self.E_total
            deltaE = current_energy - old_energy
            
            if current_energy > old_energy:
                if self.T != 0.0:
                    p = np.exp(-deltaE/(kB*self.T))
                    if random.random() < p:
                        old_energy = current_energy
                    else:
                        self.objects[index].mag = old_mag
            else:
                old_energy = current_energy
                
        if display:
            print("Monte Carlo done.\n")
        
        self.update_energies()
        self.update_caracteristics()







































































