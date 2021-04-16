# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8


import numpy as np
import matplotlib.pyplot as plt
    

def dipolar_int2D(object1, object2, angle1, angle2):
    r = np.array(object1.pos) - np.array(object2.pos)
    d = np.linalg.norm(r)
    
    if d==0:
        return 0
    
    u = r/d
    
    mag1 = [np.cos(angle1), np.sin(angle1)]
    mag2 = [np.cos(angle2), np.sin(angle2)]
    
    Ed = 3*np.dot(mag1, u)*np.dot(mag2, u)- np.dot(mag1, mag2)
    return -1e-7 * Ed / (d*1e-9)**3 * object1.M * object2.M


def dipolar_int3D(object1, object2, theta1, phi1, theta2, phi2):
    r = np.array(object1.pos) - np.array(object2.pos)
    d = np.linalg.norm(r)
    
    if d==0:
        return 0
    
    u = r/d
    
    mag1 = [np.sin(theta1)*np.cos(phi1),
           np.sin(theta1)*np.sin(phi1),
           np.cos(theta1)]
    
    mag2 = [np.sin(theta2)*np.cos(phi2),
           np.sin(theta2)*np.sin(phi2),
           np.cos(theta2)]
    
    Ed = 3*np.dot(mag1, u)*np.dot(mag2, u)- np.dot(mag1, mag2)
    return -1e-7 * Ed / (d*1e-9)**3 * object1.M * object2.M


def raw_dipolar_int3D(mag1, mag2, object1, object2):
    r = np.array(object1.pos) - np.array(object2.pos)
    d = np.linalg.norm(r)
    
    if d==0:
        return 0
    
    u = r/d
    
    Ed = 3*np.dot(mag1, u)*np.dot(mag2, u)- np.dot(mag1, mag2)
    return -1e-7 * Ed / (d*1e-9)**3 * object1.M * object2.M



def total_energy2D(angles, system, objects):
    #to minimize, in eV
    
    E=0
    
    for i in range(len(objects)):
        mag = [np.cos(angles[i]), np.sin(angles[i])]
        
        # shape
        if objects[i].ku != 0 and not objects[i].frozen and not objects[i].atomic and system.dipolar:
            theta = angles[i] - objects[i].angle * np.pi/180  
            E += objects[i].V * objects[i].ku * np.sin(theta)**2
        
        # biaxial
        if system.K1 != 0 and not objects[i].frozen:
            dots = np.dot(mag, system.ea1)**2 * np.dot(mag, system.ea2)**2
            E += objects[i].V * system.K1 * dots
            
        # uniaxial
        if objects[i].atomic and objects[i].Ku != 0:
            ku = objects[i].Ku
            if ku > 0: 
                E += ku * objects[i].V * (1-(np.dot(mag, objects[i].u_axis))**2)
            else:
                E += -ku * objects[i].V * (np.dot(mag, objects[i].u_axis))**2
        
        
        # dipolar
        if system.dipolar:
            for j in range(i):
                
                angle1 = angles[i]
                angle2 = angles[j]
                
                if objects[i].frozen:
                    mag1 = objects[i].mag
                    angle1 = np.arctan2(mag1[1], mag1[0])
                
                if objects[j].frozen:
                    mag2 = objects[j].mag
                    angle2 = np.arctan2(mag2[1], mag2[0])
                    
                E += dipolar_int2D(objects[i], objects[j], angle1, angle2)
        
        
        # Zeeman
        if system.B != 0 and not objects[i].frozen:
            m = objects[i].M*np.array([np.cos(angles[i]), np.sin(angles[i])])
            B = system.B*np.array(system.Bu)
            E += -np.dot(m, B)
        
    # exchange
    for i in range(len(system.couples)):
        angle1 = angles[system.couples[i].index1]
        angle2 = angles[system.couples[i].index2]
        
        mag1 = [np.cos(angle1), np.sin(angle1)]
        mag2 = [np.cos(angle2), np.sin(angle2)]
        
        if objects[system.couples[i].index1].frozen:
            mag1 = objects[system.couples[i].index1].mag
            
        if objects[system.couples[i].index2].frozen:
            mag2 = objects[system.couples[i].index2].mag
        
        E += -system.couples[i].J*np.dot(mag1, mag2)*1.602e-19
    
    return E/1.602e-19 #eV
    

def total_energy_min2D(angles, *args):
    # bypass the nonsense
    return total_energy2D(angles, args[0], args[0].objects)





def total_energy3D(angles, system, objects):
    #to minimize, in eV
    E=0.0
    
    for i in range(len(objects)):
        theta = angles[2*i]
        phi = angles[2*i+1]
        
        mag = [np.sin(theta)*np.cos(phi),
               np.sin(theta)*np.sin(phi),
               np.cos(theta)]
               
        
        # cubic
        if system.K1 != 0 and not objects[i].frozen:
            a = np.dot(mag, system.ea1)
            b = np.dot(mag, system.ea2)
            c = np.dot(mag, system.ea3)
            dots = (a*b)**2 + (b*c)**2 + (c*a)**2
            E += system.objects[i].V * system.K1 * dots
            
        # uniaxial
        if objects[i].atomic and objects[i].Ku != 0:
            ku = objects[i].Ku
            if ku > 0: 
                E += ku * objects[i].V * (1-(np.dot(mag, objects[i].u_axis))**2)
            else:
                E += -ku * objects[i].V *np.dot(mag, objects[i].u_axis)**2
            
        # shape
        if not objects[i].frozen and not objects[i].atomic and system.dipolar:
            phi = objects[i].angle*np.pi/180
            mag_eff = [mag[0]*np.cos(phi) + mag[1]*np.sin(phi),
                       mag[0]*np.sin(phi) - mag[1]*np.cos(phi),
                       mag[2]]
            mag_eff /= np.linalg.norm(mag_eff)
            
            demag = mag_eff[0]**2*objects[i].factors[0]
            demag += mag_eff[1]**2*objects[i].factors[1]
            demag += mag_eff[2]**2*objects[i].factors[2]
            E += system.objects[i].M0**2 * objects[i].V *(2*np.pi*1e-7)*demag
        
        # dipolar
        if system.dipolar:
            for j in range(i):
                
                theta1 = angles[2*i]
                phi1 = angles[2*i+1]
                
                theta2 = angles[2*j]
                phi2 = angles[2*j+1]
                
                if objects[i].frozen:
                    theta1 = np.arccos(objects[i].mag[2])
                    phi1 = np.arctan2(objects[i].mag[1], objects[i].mag[0])
                
                if objects[j].frozen:
                    theta2 = np.arccos(objects[j].mag[2])
                    phi2 = np.arctan2(objects[j].mag[1], objects[j].mag[0])
                
                E += dipolar_int3D(system.objects[i], system.objects[j], 
                                   theta1, phi1, theta2, phi2)
        
        
        # Zeeman
        if system.B != 0 and not objects[i].frozen:
            m = system.objects[i].M * np.array(mag)
            B = system.B*np.array(system.Bu)
            E += -np.dot(m, B)
        
    # exchange
    for i in range(len(system.couples)):
        i1 = system.couples[i].index1
        i2 = system.couples[i].index2
        
        theta1 = angles[2*i1]
        phi1 = angles[2*i1+1]
        
        mag1 = [np.sin(theta1)*np.cos(phi1),
               np.sin(theta1)*np.sin(phi1),
               np.cos(theta1)]
        
        theta2 = angles[2*i2]
        phi2 = angles[2*i2+1]
        
        mag2 = [np.sin(theta2)*np.cos(phi2),
               np.sin(theta2)*np.sin(phi2),
               np.cos(theta2)]
        
        if objects[system.couples[i].index1].frozen:
            mag1 = objects[system.couples[i].index1].mag
            
        if objects[system.couples[i].index2].frozen:
            mag2 = objects[system.couples[i].index2].mag
        
        E += -system.couples[i].J*np.dot(mag1, mag2)*1.602e-19
    
    return E/1.602e-19 #eV
    

def total_energy_min3D(mags, *args):
    return total_energy3D(mags, args[0], args[0].objects)


def field_energy(mag, system, obj):
    #to minimize, in eV
    E=0.0
    
    # cubic
    if system.K1 != 0:
        a = np.dot(mag, system.ea1)
        b = np.dot(mag, system.ea2)
        c = np.dot(mag, system.ea3)
        dots = (a*b)**2 + (b*c)**2 + (c*a)**2
        E += obj.V * system.K1 * dots
    
    # uniaxial
    if obj.atomic and obj.Ku != 0:
        if obj.Ku > 0: 
            E += obj.Ku * obj.V * (1-(np.dot(mag, obj.u_axis))**2)
        else:
            E += -obj.Ku * obj.V * (np.dot(mag, obj.u_axis))**2
        
    # shape
    if not obj.atomic and system.dipolar:
        phi = obj.angle*np.pi/180
        mag_eff = [mag[0]*np.cos(phi) + mag[1]*np.sin(phi),
                   mag[0]*np.sin(phi) - mag[1]*np.cos(phi),
                   mag[2]]
        mag_eff /= np.linalg.norm(mag_eff)
        
        demag = mag_eff[0]**2*obj.factors[0]
        demag += mag_eff[1]**2*obj.factors[1]
        demag += mag_eff[2]**2*obj.factors[2]
        E += obj.M0**2 * obj.V *(2*np.pi*1e-7)*demag
    
    # dipolar
    '''
    if system.dipolar:
        for j in range(len(system.objects)):
            if obj.pos == system.objects[j].pos:
                continue
            
            E += raw_dipolar_int3D(mag, system.objects[j].mag,
                               obj, system.objects[j]) / 2.0
    '''
    
    # Zeeman
    if system.B != 0:
        B = system.B*np.array(system.Bu)
        E += - obj.M * np.dot(mag, B)
    
    # exchange
    for k in range(len(system.couples)):
        i1 = system.couples[k].index1
        if mag == system.objects[i1].mag:
            
            i2 = system.couples[k].index2
            mag2 = system.objects[i2].mag
            
            E += -system.couples[k].J*np.dot(mag, mag2)*1.602e-19

    return E


def field_energy_bypass(mag, *args):
    return field_energy(mag, args[0], args[1])











