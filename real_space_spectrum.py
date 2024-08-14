# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 17:53:41 2024

@author: harry
"""

from magnetic_interface_class_file import *


if __name__ == "__main__":
    
    Vm =    1
    phi =   0
    theta = 0

    params = {
        'Nx'            : 30,
        'Ny'            : 30,
        't'             : 1,
        'Delta'         : 0.1,
        'mu'            : -3.6,
        'km'            : 0.65,
        'Nev'           : 100,
        'bc'            : 'periodic',
        'impurities'    : []
    }

    
    
    sc = SC_localiser(**params)
    
    Vm_values=np.linspace(0,6,51)
    
    spectrum=np.zeros((params["Nev"],len(Vm_values)))
    
    for Vm_indx,Vm in enumerate(tqdm(Vm_values)):
        params["Vm"]=Vm
        impurities=circle_spiral(params, params["Nx"]//4, 50, "radius", Vm, phi, theta)
        params["impurities"]=impurities
        sc.Create_Hamiltonian()
        sc.Solve_Minimal()
        spectrum[:,Vm_indx]=sc.En
        
    plt.figure()
    for i in range(params["Nev"]):
        plt.plot(Vm_values,spectrum[i,:],"kx")

        
        
        