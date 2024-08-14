# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:13:59 2024

@author: Harry MullineauxSanders
"""

from Circle_function import *
import matplotlib.pyplot as plt
Vm =    1
phi =   0
theta = 0

params = {
    'Nx'            : 200,
    'Ny'            : 200,
    't'             : 1,
    'Delta'         : 0.1,
    'mu'            : -3.6,
    'km'            : 0.65,
    'Nev'           : 2,
    'bc'            : 'periodic',
    'impurities'    : []
}

impurities=circle_spiral(params, 100, 50, "radius", Vm, phi, theta)

plt.figure()

for impurity in impurities:
    plt.plot(impurity[0],impurity[1],"x")