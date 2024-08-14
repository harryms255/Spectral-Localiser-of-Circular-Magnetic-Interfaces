# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:07:09 2024

@author: Harry MullineauxSanders
"""

from time import time
import numpy as np
import scipy.sparse.linalg as sl
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.sparse import dok_matrix
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from Circle_function import *
plt.close("all")
plt.rcParams.update({'font.size': 30})
num_colors = 3
cmap = colors.ListedColormap(['blue', 'yellow'])


def line_spiral(params, l, Vm, phi, theta):
    
    '''
    This functions creates a line of impurites. The spin orientation is helical.
    
    Parameters
    ----------
    
    params : dict
        dictionary containing system profile
    l : int
        length of chain
    Vm : float
        magnetic scattering strength of impurities
    phi : float
        azimuthal angle
    theta : float
        polar angle
    
    Returns
    -------
        array of characteristics for each impurity
    
    '''
    # creates a line of impurities
    impurities = []

    Nx2 = params['Nx']//2
    Ny2 = params['Ny']//2

    # optimal pitch shift
    dphi = 2*np.arccos((-2*params['t'] - params['mu'])/(2*params['t']))

    for idx in range(l):
        xn = (Nx2 - l//2) + idx
        impurities.append([xn, Ny2, Vm, idx*dphi, theta, 0.])
    return impurities


def conjT(U):
    return np.conjugate(U).T


class SC_localiser(object):

    '''
    Object Set-up
    '''
    def __init__(self, *args, **kwargs):
        self._parameters = {
            'pauli_x'         : sp.csc_matrix(np.array(([0,1],[1,0]), dtype=complex)),
            'pauli_y'         : sp.csc_matrix(np.array(([0,-1j],[1j,0]), dtype=complex)),
            'pauli_z'         : sp.csc_matrix(np.array(([1,0],[0,-1]), dtype=complex))
        }
        
        self._parameters.update(kwargs)
        
        for key in list(self._parameters.keys()):
            setattr(self, key, self._parameters[key])

        self.tic = time()
    
        self.Update_Parameters(**kwargs)
    
    def __enter__(self):
        return self
    
    def __exit__(self,*err):
        pass
    
    def __del__(self):
        pass
    
    def Update_Parameters(self, **kwargs):
        
        self._parameters.update(kwargs)
        
        for key in list(self._parameters.keys()):
            setattr(self, key, self._parameters[key])
        
    '''
    Functions
    '''
    
    def Load_Hamiltonian(self):
        path = '/Users/Christina/Desktop/Project/'
        self.H = (np.load(path + 'H_150_50.npy', allow_pickle=True)).tolist()
        
    # -------------------------------------------------------------------------

    def Solve_All(self, create = True):

        En0, U0 = sl.eigsh(self.H, 2*self.Nx*self.Ny, which='SM')
        En1, U1 = sl.eigsh(self.H, 2*self.Nx*self.Ny, sigma=0, which='SM')
 
        self.En = np.concatenate([En0,En1])
        self.U = np.concatenate([U0,U1],axis=1)
        
    # -------------------------------------------------------------------------

    def Solve_Minimal(self, create = True):
        
        self.En, self.U = sl.eigsh(self.H, self.Nev, sigma=0, which ='SM')

    # -------------------------------------------------------------------------
    def x_operator(self):
        self.x_hat = sp.csc_matrix(np.diag(np.tile(np.repeat(np.arange(self.Nx), 4)%self.Nx, self.Ny)))
    
    # -------------------------------------------------------------------------
    
    def y_operator(self):
        self.y_hat = sp.csc_matrix(np.diag(np.repeat(np.repeat(np.arange(self.Ny),4)%self.Ny, self.Nx)))
    
        
    
    #--------------------------------------------------------------------------
    def angle_operator(self):
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
            
        self.theta_hat
    
    
    # -------------------------------------------------------------------------
    
    def Lattice_size(self):
        return self.Nx, self.Ny
    
    # -------------------------------------------------------------------------
    
    def spectral_localiser(self, x, y,E):
        #Function that calculates the Spectral Localiser Operator for a given Hamiltonian at a given position and Energy
        #Form of the localiser taken from PHYSICAL REVIEW B 106, 064109 (2022)
        
        k=0.01
        spectral_localiser=k*(self.X-x*self.iden_X+self.Y-y*self.iden_Y)+self.Z-E*self.iden_Z
        
        return spectral_localiser

    # -------------------------------------------------------------------------

    
    def kron(self):
        
        """
        In the spectral localiser we make use of tensor products of large matricies
        so it is efficent to do this only once at the start
        """
        self.X=sp.kron(self.pauli_x,self.x_hat)
        self.Y=sp.kron(self.pauli_y,self.y_hat)
        self.Z=sp.kron(self.pauli_z,self.H)
        
        self.iden_X=sp.kron(self.pauli_x,sp.identity(4*self.Nx*self.Ny))
        self.iden_Y=sp.kron(self.pauli_y,sp.identity(4*self.Nx*self.Ny))
        self.iden_Z=sp.kron(self.pauli_z,sp.identity(4*self.Nx*self.Ny))
    
    # -------------------------------------------------------------------------
    
    def localiser_gap(self, x, y,E):
        #Function that calculates the gap in the spectrum of the spectral localiser
        #Diagonalises the full matrix and returns the smallest eigenvalue
        #This seems to be quicker than calculating the lowest eigenvalue and eigenvector directly via spare matrix operations

        eigenvalues = sl.eigsh(self.spectral_localiser(x, y,E), self.Nev, sigma=0, which ='SM', return_eigenvectors=False)
        
        #eigenvalues=np.linalg.eigvalsh(self.spectral_localiser(x, y, E).A)
        return np.min(abs(eigenvalues))
    
    # -------------------------------------------------------------------------
    
    def class_D_invariant(self, x):
        #TERRY A LORING, "K-theory and pseudospectra for topological insulators", ANNALS OF PHYSICS, 356, 2016, 383-416
        
        y = self.Ny // 2
        
        C = self.x_hat - x * sp.identity(4 * self.Nx * self.Ny) + self.y_hat - y * sp.identity(4 * self.Nx * self.Ny) + 1j * self.H

        sgn, log = np.linalg.slogdet(C)

        return np.real(sgn)

    # -------------------------------------------------------------------------

    def Local_Chern_Marker(self, ph = 0, n = None):
        idx = self.En < 0
        Up  = U[:, idx]
        
        if n != None:
            Up = Up[:, -min(len(idx), n):]
    
        P = Up @ conjT(Up)
        XY = -2j*np.pi*(P @ self.x_hat @ P @ self.y_hat @ P)
        YX = -2j*np.pi*(P @ self.y_hat @ P @ self.x_hat @ P)
        #Z = XY + conjT(XY)
        Z=XY-YX

        return np.diag(Z)[ph::4].reshape(self.Ny, self.Nx).real

    # -------------------------------------------------------------------------
    
    def Create_Hamiltonian(self):
        
        Nx, Ny = self.Lattice_size()
        
        H=dok_matrix((4*Nx*Ny,4*Nx*Ny),dtype=complex)
        
        #over all y sites
        for ny in range(self.Ny):

            # loop over all x sites
            for nx in range(self.Nx):

                # lattice index
                n = nx + ny*self.Nx

                # index to matrix entries (4 * larger)
                n4 = 4*n

                # ------------------------------------------------------------

                damping_bond_x = 1.0
                damping_bond_y = 1.0

                # set hopping amplitudes with the damping
                tx = self.t * damping_bond_x
                ty = self.t * damping_bond_y

                # ------------------------------------------------------------

                # on-site: subtract/add chemical potential
                #
                # The Hamiltonian has the following structure
                #
                # H = [ (band 1)  (gaps)   ]
                #     [ (gaps)^+ -(band 2) ]
                #
                # where (band 1) and (band 2) are the 2x2 Hamiltonians for the normal
                # contributions to the Hamitonian, written in the particle sector
                # (~ c^+ c) as "band 1" and for the hole sector (~ c c^+) as "band 2"
                # which provides for the latter all the - signs for the energies.
                #
                # The arrangement of basis vectors is
                #
                #    [ c_{1,up} c_{1,dw} c_{2,up}^+ c_{2,dw}^+ ]
                #
                # there's no minus sign in this vector (which is occasionally used in
                # the literature).

                # diagonal parts of Hamiltonian, minus signs for band 2
                # because it's in the hole representation
                H[n4,  n4  ] = -self.mu
                H[n4+1,n4+1] = -self.mu
                H[n4+2,n4+2] =  self.mu
                H[n4+3,n4+3] =  self.mu


                # ------------------------------------------------------------
                # on-site gap pairing

                # For the structure mentioned already above
                #
                # H = [ (band 1)   (gaps)   ]
                #     [ (gaps)^+  -(band 2) ]
                #
                # we have the gaps
                #
                # gaps = [ Delta(up,up) Delta(up,dw) ]
                #        [ Delta(dw,up) Delta(dw,dw) ]
                #
                # corresponding to the entries in the Hamiltonian:
                #
                #        [ (n4,  n4+2) (n4,  n4+3) ]
                #        [ (n4+1,n4+2) (n4+1,n4+3) ]
                #
                # only the singlet sector will be set:
                #
                # Delta(up,dw) = - Delta(dw,up) = Delta

                # upper right block "(gaps)"
                H[n4,  n4+3] =  self.Delta
                H[n4+1,n4+2] = -self.Delta

                # fill remaining entries., lower left block, "(gaps)^+", with h.c.
                H[n4+2,n4+1] = H[n4+1,n4+2].conj()
                H[n4+3,n4  ] = H[n4,  n4+3].conj()


                # ------------------------------------------------------------
                # hopping between neighbouring x sites:
                # - choose neighbour mx = nx + 1 (modulo Nx for periodic bc)
                # - m = mx + my*Ny
                # - entries in block 4x4 matrix indexed by (n,m):
                #   [ ty       0       0       0   ]
                #   [ 0        ty      0       0   ]
                #   [ 0        0      -ty      0   ]
                #   [ 0        0       0      -ty  ]

                # neighbour along +x (the modulo % ensures periodic boundary cond)
                # and -x neighbour by h.c.
                if (self.bc == "periodic") or (nx+1 < self.Nx):
                    mx = (nx+1) % self.Nx
                    my = ny
                    m  = mx + my*self.Nx
                    m4 = 4*m

                    # --------------------------------------------------------

                    # hopping: only diagonal entries in the 4x4 block as without
                    # spin-orbit interaction the hopping is spin preserving
                    H[n4,  m4  ] =  tx
                    H[n4+1,m4+1] =  tx
                    H[n4+2,m4+2] = -tx
                    H[n4+3,m4+3] = -tx

                    # block matrix for (ny+1,ny) by h.c. (the -y neighbour for ny+1)
                    H[m4,  n4  ] = H[n4,  m4  ].conj()
                    H[m4+1,n4+1] = H[n4+1,m4+1].conj()
                    H[m4+2,n4+2] = H[n4+2,m4+2].conj()
                    H[m4+3,n4+3] = H[n4+3,m4+3].conj()


                # ------------------------------------------------------------
                # hopping between neighbouring y sites:
                # - choose neighbour my = ny + 1 (modulo Ny for periodic bc)
                # - entries in block 4x4 matrix indexed by (ny,my):
                #   [ ty       0       0       0  ]
                #   [ 0        ty      0       0  ]
                #   [ 0        0      -ty      0  ]
                #   [ 0        0       0      -ty ]

                # neighbour along +y (the modulo % ensures periodic boundary cond)
                # and -y neighbour by h.c.
                if (self.bc == "periodic") or (ny+1 < self.Ny):
                    mx = nx
                    my = (ny+1) % self.Ny
                    m  = mx + my*self.Nx
                    m4 = 4*m

                    # --------------------------------------------------------

                    # hopping: only diagonal entries in this 4x4 block
                    H[n4,  m4  ] =  ty
                    H[n4+1,m4+1] =  ty
                    H[n4+2,m4+2] = -ty
                    H[n4+3,m4+3] = -ty

                    # block matrix for (ny+1,ny) by h.c. (the -y neighbour for ny+1)
                    H[m4,  n4  ] = H[n4,  m4  ].conj()
                    H[m4+1,n4+1] = H[n4+1,m4+1].conj()
                    H[m4+2,n4+2] = H[n4+2,m4+2].conj()
                    H[m4+3,n4+3] = H[n4+3,m4+3].conj()


        # add impurities
        for imp in self.impurities:
            # each entry in `impurities` is of the form
            #  ( nx, ny, Vm, Vm_phi, Vm_theta, Vs )
            x,y,Vm,Vm_phi,Vm_theta,Vs = imp

            n = x + y*self.Nx
            n4 = 4*n

            Vmx = Vm * np.cos(Vm_theta) * np.cos(Vm_phi)
            Vmy = Vm * np.cos(Vm_theta) * np.sin(Vm_phi)
            Vmz = Vm * np.sin(Vm_theta)

            # for the Hamiltonian we nees Vplus,Vminus
            Vmp = Vmx + 1.j * Vmy
            #Vmm = Vmx - 1.j * Vmy


            # update Hamiltonian (make sure to use += and not
            # = here to avoid overwriting what's been added above)

            # band 1: particle representation
            # up up: + Vm^z and scalar potential
            H[n4,  n4  ] +=  Vmz + Vs
            # up dw: Vm^+
            H[n4,  n4+1]  =  Vmp.conj()
            # dw up: Vm^- [ = (Vm+^)^* ]
            H[n4+1,n4  ]  =  Vmp
            # dw dw: - Vm^z and scalar potential
            H[n4+1,n4+1] += -Vmz + Vs

            # band 2: hole representation: extra minus signs to all
            # quantities
            H[n4+2,n4+2] += -Vmz - Vs
            H[n4+2,n4+3]  = -Vmp
            H[n4+3,n4+2]  = -Vmp.conj()
            H[n4+3,n4+3] +=  Vmz - Vs

        # if sparse convert to csc format
        self.H = H.tocsc()




