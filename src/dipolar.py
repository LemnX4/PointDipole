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


def dipolar_field(position, obj):
    r = np.array(position) - np.array(obj.pos)
    d = np.linalg.norm(r)
    
    u = r/d
    
    field = 3*np.dot(obj.mag, u)*u - obj.mag
    return -1e-7/(d*1e-9)**3 * field * obj.M




















