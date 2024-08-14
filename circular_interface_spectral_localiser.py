# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:11:51 2024

@author: Harry MullineauxSanders
"""

from magnetic_interface_class_file import *

if __name__ == "__main__":
    
    Vm =    1
    phi =   0
    theta = 0

    params = {
        'Nx'            : 20,
        'Ny'            : 20,
        't'             : 1,
        'Delta'         : 0.1,
        'mu'            : -3.6,
        'km'            : 0.65,
        'Nev'           : 2,
        'bc'            : 'periodic',
        'impurities'    : []
    }

    impurities=circle_spiral(params, params["Nx"]//4, 50, "radius", Vm, phi, theta)
    
    sc = SC_localiser(**params)
    
    sc.x_operator()
    sc.y_operator()
    sc.Create_Hamiltonian()
    sc.kron()
    
    
    x_values=np.linspace(0,params["Nx"]-1,params["Nx"])
    y_values=np.linspace(0,params["Ny"]-1,params["Ny"])
    
    
    localiser_gap_values=np.zeros((len(y_values),len(x_values)))
    
    for x_indx,x in enumerate(tqdm(x_values)):
        for y_indx,y in enumerate(tqdm(y_values)):
            localiser_gap_values[y_indx,x_indx]=sc.localiser_gap(x, y,0)
    
    fig,ax=plt.subplots()
    sns.heatmap(localiser_gap_values,cmap="plasma",vmin=0,ax=ax)
    ax.invert_yaxis()
    
    
    