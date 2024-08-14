# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 15:31:18 2024

@author: Harry MullineauxSanders
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.close("all")

def theta(x,y,x0,y0):
    if (x-x0)>0:
        theta=np.arctan((y-y0)/(x-x0))
    elif (y-y0)>=0:
        if (x-x0)==0:
            theta=np.pi/2
        else:
            theta=np.pi+np.arctan((y-y0)/(x-x0))
    elif (y-y0)<0:
        if (x-x0)==0:
            theta=-np.pi/2
        else:
            theta=-np.pi+np.arctan((y-y0)/(x-x0))
    elif x==0:
        theta=0
        
    return theta
    
   

Nx=101
Ny=101
x0=Nx//2
y0=Ny//2
x_values=np.linspace(0,Nx-1,1001)
y_values=np.linspace(0,Ny-1,1001)

theta_values=np.zeros((len(y_values),len(x_values)))

for y_indx,y in enumerate(y_values):
    for x_indx,x in enumerate(x_values):
        theta_values[y_indx,x_indx]=theta(x, y, x0, y0)

plt.figure()
sns.heatmap(theta_values/np.pi,cmap="viridis",vmin=-1,vmax=1)
plt.gca().invert_yaxis()