# acoustic pressure modelisation
# based on resulstut from literature_comp_PETRA_acoustic_pressure.py
# use equation (7) with (13) from [1] 
# called P_bis(r,t) in the code literature_comp_PETRA_acoustic_pressure.py 
##### REF :
# [1] Petra, N., Zweck, J., Kosterev, A. A., Minkoff, S. E., & Thomazy, D. (2009). 
# Theoretical analysis of a quartz-enhanced photoacoustic spectroscopy sensor. 
# Applied Physics B: Lasers and Optics, 94(4), 673–680. 
# https://doi.org/10.1007/s00340-009-3379-1
# /!\ see also pdf and matlab file from her student
##### VERSION
# 2020-04-17 - first version

from scipy.special import j0, y0 # zeroth-order Bessel functions of the first and second kinds
from scipy.integrate import quad # integration 
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
# get current working directory
path = os.path.abspath(os.getcwd())
index = path[::-1].find('/')
path = path[:len(path)-index-1]
#print('current working directory : ',path)
# set the path of the module I want
sys.path.append(path) 

from fluid_material import fluid_material as fm

class AcousticPressure():
    # define a class wich can return the acoustic pressure
    def __init__(self,fluid):
        self.fluid = fluid
        # frequency
        self.freq = 32761 #[Hz]        
        # laser beam width
        self.sigma = 0.05e-3 #[m]
        # laser power
        self.Pl = 61.7e-3 #[W]
        # effecrive losses
        self.k_eff = 1.31e-2 #[m-1]

    def P(self,r,t):
        # equation (7) with (13) from [1]
        # laser beam width
        sigma = self.sigma  #[m]
        omega = 2*np.pi*self.freq #[1/s]
        # air heat capacity ratio
        Gamma = self.fluid.Gamma
        # amplitude of the acoustic source
        W = -1*(Gamma-1)*self.k_eff/2*self.Pl/(2*np.pi*self.sigma**2)*omega
        # print('Amplitude of the acoustic source : %g N/m2/s2'%W)
        # speed of sound
        c = self.fluid.speed_of_sound
        def f1_minus_i_f2(r):
            # equation (13) from [1]    
            k = omega/c
            A = -np.pi*W/2./c**2/k**2
            def h(u):
                return u*j0(u)*np.exp(-u**2/2./k**2/sigma**2)
            int_h, _ = quad(h, 0., np.inf)
            return A*int_h*(y0(k*r)+1j*j0(k*r))
        pr = f1_minus_i_f2(r)
        return pr*np.exp(1j*omega*t)

if __name__ == "__main__":
    Air = fm.Fluid(name='air')
    Air.info()
    AcPress = AcousticPressure(Air)

    r = np.linspace(1e-9,2e-3,1000)
    P = np.asarray([AcPress.P(i,0) for i in r])

    # graph
    gridsize = (1, 2)
    fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
    fig.suptitle('PETRA Acoustic pressure', fontsize=20)

    ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
    ax1.plot(r*1e3, np.abs(P)*1e3,linewidth=2,label='equation (7) with (13)')
    ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.set_xlabel('Radial distance from laser beam (mm)', labelpad=5, fontsize=14)
    ax1.set_ylabel('Pressure amplitude (mPa)', labelpad=5, fontsize=14)
    # ax1.set_title('title',fontsize=16)
    ax1.set_xlim([0,2])
    ax1.set_ylim([0,0.2])
    ax1.legend(loc=3)

    ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
    ax2.plot(r*1e3, np.angle(P,deg=True),linewidth=2, label='equation (7) with (13)')
    ax2.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.set_xlabel('Radial distance from laser beam (mm)', labelpad=5, fontsize=14)
    ax2.set_ylabel('Pressure phase (degrees)', labelpad=5, fontsize=14)
    # ax2.set_title('title',fontsize=16)
    # ax2.set_xlim([0,2])
    # ax2.set_ylim([0,120])
    ax2.legend(loc=0)
    plt.show()
