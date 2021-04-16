# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8


import numpy as np
from scipy.optimize import approx_fprime

from energy import field_energy_bypass
from dipolar import dipolar_field

def Bth(gamma, alpha, M, T):
    eta = [np.random.normal(0, 1), np.random.normal(0, 1), np.random.normal(0, 1)]
    kB = 1.380649e-23
    Bt = np.array(eta)*np.sqrt(2*alpha*kB*T/(M*gamma))
    return Bt


def LLGfunc(M, B, Ms, gamma, alpha):
    # dM = LLGfunc * dt
    return -gamma*(np.cross(M, B) + alpha/Ms * np.cross(M, np.cross(M, B)))


def LLG_total(t, mags, system, gamma, alpha, dt, epsilon):
    new_mags = []
        
    for i in range(len(system.objects)):
        mag = np.array([mags[3*i], mags[3*i+1], mags[3*i+2]])
        
        epsilons = [epsilon, epsilon, epsilon]
        B = -approx_fprime(mag, field_energy_bypass, epsilons, 
                          system, system.objects[i])/system.objects[i].M
        
        if system.dipolar:
            for j in range(len(system.objects)):
                if i==j:
                    continue
                B += dipolar_field(system.objects[i].pos, system.objects[j])
                           
        if system.T != 0:
            B += system.objects[i].thermic_field[int(t/dt)]/np.sqrt(dt)
        
        new_mag = LLGfunc(mag*system.objects[i].M, B, system.Ms, gamma, alpha)/system.objects[i].M
        
        new_mags.append(new_mag[0])
        new_mags.append(new_mag[1])
        new_mags.append(new_mag[2])
    
    return new_mags

def LLG(t, mags, *args):
    system = args[0]
    gamma = args[1]
    alpha = args[2]
    dt = args[3]
    epsilon = args[4]
    
    return LLG_total(t, mags, system, gamma, alpha, dt, epsilon)













