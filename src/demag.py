# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8


import numpy as np
import scipy.integrate as integrate
import scipy.special as special

import matplotlib.pyplot as plt



def g(beta, theta):
    return 1-(1-beta**2)*np.cos(theta)**2


def F(x):
    return special.hyp2f1(-0.5, 0.5, 2, x)


def Nz(beta, tau):
    
    A = 1 + 8*beta/(3*np.pi**2*tau)*special.ellipk(1-beta**2)
    
    integral = integrate.quad(lambda theta: F(-beta**2/(tau**2*g(beta, theta))), 0, np.pi/2)
    
    return A - 2/np.pi * integral[0]


def Nx(beta, tau):
    
    A = 8*beta/(3*np.pi**2 *tau*(1-beta**2)) * (beta**2 *special.ellipk(1-beta**2) - special.ellipe(1-beta**2))
    
    integral = integrate.quad(lambda theta: np.cos(theta)**2/g(beta, theta)*F(-beta**2/(tau**2*g(beta, theta))), 0, np.pi/2)
    
    return A + 2*beta**2/np.pi * integral[0]


def demagnetization(a, b, h, show=False):
    nz = Nz(b/a, h/a)
    
    if a==b:
        nx = (1 - nz)/2
        ny = nx
    else:
        nx = Nx(b/a, h/a)
        ny = 1 - nx - nz
    
    factors = [nx, ny, nz]
    
    if show:
        print("\nNx = " + str(nx))
        print("Ny = " + str(ny))
        print("Nz = " + str(nz))
        print("\nMinimun energy axis : " + (["Ox", "Oy", "Oz"])[factors.index(min(factors))] + "\n")
    
    return factors


def ku(a, b, h, M0):
    factors = demagnetization(a, b, h)
    
    return M0**2 *(2*np.pi*1e-7)*(1 - 2*factors[0] - factors[2])



if __name__ == "__main__":
    
    a = 4e-9
    b = 4e-9
    h = 3.5e-9
    M0 = 1720e3
    
    
    
    
    sizes = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]
    
    energies = []
    
    for i in range(len(sizes)):
        a = sizes[i]*1e-9
        b = a/1.2
        h = a/2.0
        M0 = 1720e3
        
        V = np.pi * (a/2.0)*(b/2.0) * h
    
        factors = demagnetization(a, b, h)
        
        E0 = 0.5 * (4*np.pi*1e-7) * M0**2 * V * factors[0]
    
    
        print("E0= " + str(E0) + " J")
        
        energies.append(E0)
    
    
    plt.plot(sizes, energies)
    plt.xlabel("Size (nm)")
    plt.ylabel("E0 (J)")
    plt.title("Internal dipolar energy for different sizes (VAR=1.2)")
    plt.savefig("internal_energy.png")#, dpi=200)
    plt.clf()
    
    
    '''
    ars = [1.00001, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09,
          1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19,
          1.2]
    
    constants =  []
    
    for i in range(len(ars)):
        a = 4e-9
        b = a/ars[i]
        h = 2e-9
        M0 = 1720e3
        constants.append(Ku(a, b, h, M0)/1000)
    
    
    plt.plot(ars, constants)
    plt.xlabel("Aspect ratio")
    plt.ylabel("Uniaxial constant (kJ/m3)")
    plt.title("Theoritical uniaxial constant for different aspect ratios.")
    plt.savefig("theory_constants_horizontal.png")#, dpi=200)
    plt.clf()
    
    
    
    ars = [1.00001, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]
    
    constants =  []
    
    for i in range(len(ars)):
        a = 4e-9
        b = a/1.1
        h = 4e-9/ars[i]
        M0 = 1720e3
        constants.append(Ku(a, b, h, M0)/1000)
    
    
    plt.plot(ars, constants)
    plt.xlabel("Vertical aspect ratio")
    plt.ylabel("Uniaxial constant (kJ/m3)")
    plt.title("Theoritical uniaxial constant for different aspect ratios.")
    plt.savefig("theory_constants_vertical.png")#, dpi=200)
    plt.clf()
    
    '''



















