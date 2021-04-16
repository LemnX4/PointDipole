# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8


import numpy as np
import random

kB = 1.380649e-23


def normalize(p):
    part = 0.0
    for i in range(len(p)):
        part += p[i]
    for i in range(len(p)):
        p[i] /= part

def mean(l, p):
    s = 0.0
    for i in range(len(l)):
        s += p[i]*l[i]
    return s
    
def squareMean(l, p):
    s = 0.0
    for i in range(len(l)):
        s += p[i]*l[i]**2
    return s
    
def susceptibility(m, p, T):
    return (squareMean(m, p) - mean(m, p)**2)/(kB*T)

def specificHeat(e, p, T):
    return (squareMean(e, p) - mean(e, p)**2)/(kB*T**2)

def thermal_at(system, N, NMC, T):
    p = []
    e = []
    m = []
    old_T = system.T
    for n in range(N):
        system.T = T
        system.MonteCarlo(NMC, display=False)
        
        pro = np.exp(-system.E_total/(kB*system.T))
        
        p.append(pro)
        e.append(system.E_total)
        m.append(system.M/system.Ms)
    
    system.T = old_T
    normalize(p)
    return susceptibility(np.array(m)*system.Ms, p, T), specificHeat(e, p, T), mean(m, p)
    
def thermal(system, N, Nsample, T1, T2, filename="", display=True):
    Ts = np.linspace(0, N, N)*(T2-T1)/N + T1
    xi = []
    cv = []
    ms = []
    
    if display:
        print("\nDoing thermal measurements...")
    
    for T in Ts:
        #system.relax(display=False)
        system.randomize_magnetizations()
        x, c, m = thermal_at(system, Nsample, 20+int(T/2.0), T)
        xi.append(x)
        cv.append(c)
        ms.append(m)
    
    if display:
        print("Thermal measurements done.\n")
    
    if filename != "":
        file = open(filename, "w")
        if display:
            print('\nWriting thermal measurement data in the "{}" file...'.format(filename))
        
        file.write("T\tXi\tCv\tM\n")
        for i in range(len(Ts)):
            file.write("{}\t{}\t{}\t{}\n".format(Ts[i], xi[i], cv[i], ms[i]))
        
        if display:
            print("Done writing.\n")
        file.close()
    
    
    return Ts, xi, cv, ms








































