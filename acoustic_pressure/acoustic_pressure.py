'''
acoustic pressure modelisation
based on resulstut from literature_comp_PETRA_acoustic_pressure.py
use equation (7) with (13) from [1] 
called P_bis(r,t) in the code literature_comp_PETRA_acoustic_pressure.py 
#### REF :
[1] Petra, N., Zweck, J., Kosterev, A. A., Minkoff, S. E., & Thomazy, D. (2009). 
Theoretical analysis of a quartz-enhanced photoacoustic spectroscopy sensor. 
Applied Physics B: Lasers and Optics, 94(4), 673–680. 
https://doi.org/10.1007/s00340-009-3379-1
see also pdf and matlab file from her student
#### VERSION
2020-04-17 - first version
2020-04-19 - take into account that the laser beam does not have a constant 
            width along the optical axis.
'''

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
    def __init__(self,fluid,dico):
        self.fluid = fluid
        # modualtion frequency
        self.freq_mod = dico['freq_mod'] #[Hz]        
        # laser power
        self.Pl = dico['laser_power'] #[W]
        # laser wavelength
        self.wavelength = dico['wavelength'] # [m]
        # laser waist (minimum value possible)
        # if we compare with [1] sigma = waist/2
        if dico['waist'] is None:
            self.waist = 2*self.wavelength/np.pi
        else:
            self.waist = dico['waist']
        # laser Rayleigh length
        if dico['Rayleigh_length'] is None:
            self.Rayleigh = np.pi*self.waist**2/self.wavelength
        else:
            self.Rayleigh = dico['Rayleigh_length'] 
        # effecrive losses
        self.alpha_eff = dico['alpha_eff'] #[m-1]

    def Change_waist(self,value):
        # change the value for the waist
        # and update Rayleigth length
        self.waist = value
        # laser Rayleigh length
        self.Rayleigh = np.pi*self.waist**2/self.wavelength 

    def P(self,x,y,z,t):
        '''
        - code of equation (7) with (13) from [1]
        - the gaussian profil of the laser if different 
        take into account that the laser beam does not have a 
        - the optical axis is along the x-axis
        - the laser is above the plan xy and its heigh is zL
        '''
        # radial distance
        r = np.sqrt(y**2+z**2)
        # laser beam width in x
        wL0 = self.waist
        xR = self.Rayleigh
        wL = wL0*np.sqrt(1+x**2/xR**2)
        # laser modulation
        omega = 2*np.pi*self.freq_mod #[1/s]
        # laser power
        Pl = self.Pl
        # air heat capacity ratio
        Gamma = self.fluid.Gamma
        # effecrive losses
        k_eff = self.alpha_eff
        # amplitude of the acoustic source
        W = -1*(Gamma-1)*k_eff/2*Pl*2./(np.pi*wL**2)*omega
        # print('Amplitude of the acoustic source : %g N/m2/s2'%W)
        # speed of sound
        c = self.fluid.speed_of_sound
        def f1_minus_i_f2(r):
            # equation (13) from [1]    
            k = omega/c
            A = -np.pi*W/2./c**2/k**2
            def h(u):
                return u*j0(u)*np.exp(-2.*u**2/k**2/wL**2)
            int_h, _ = quad(h, 0., np.inf)
            return A*int_h*(y0(k*r)+1j*j0(k*r))
        pr = f1_minus_i_f2(r)
        return pr*np.exp(1j*omega*t)

    def info(self):
        # print important atribut 
        print('LASER :')
        print('wavelength       = %f um'%(self.wavelength*1e6))
        print('Rayleigh length  = %f um'%(self.Rayleigh*1e6))
        print('waist            = %f um'%(self.waist*1e6))
        print('modulation       = %f kHz'%(self.freq_mod*1e-3))
        print('GAS   :')
        print("losses           = %f cm-1"%(self.alpha_eff*1e-2))


if __name__ == "__main__":
    Air = fm.Fluid(name='air')
    Air.info()
    AcPress = AcousticPressure(Air)
    AcPress.waist = 100e-6 # [m]
    AcPress.info()

    Pressure = []
    tab_z = np.linspace(1e-9,2e-3,500)
    for i in tab_z:
        temp = AcPress.P(0,0,i,0)
        Pressure.append(temp)       
    Pressure = np.array(Pressure)

    # graph
    gridsize = (1, 2)
    fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
    fig.suptitle('PETRA Acoustic pressure', fontsize=20)

    ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
    ax1.plot(tab_z*1e3, np.abs(Pressure)*1e3,linewidth=2,label='equation (7) with (13)')
    ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.set_xlabel('Radial distance from laser beam (mm)', labelpad=5, fontsize=14)
    ax1.set_ylabel('Pressure amplitude (mPa)', labelpad=5, fontsize=14)
    # ax1.set_title('title',fontsize=16)
    ax1.set_xlim([0,2])
    ax1.set_ylim([0,0.2])
    ax1.legend(loc=3)

    ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
    ax2.plot(tab_z*1e3, np.angle(Pressure,deg=True),linewidth=2, label='equation (7) with (13)')
    ax2.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.set_xlabel('Radial distance from laser beam (mm)', labelpad=5, fontsize=14)
    ax2.set_ylabel('Pressure phase (degrees)', labelpad=5, fontsize=14)
    # ax2.set_title('title',fontsize=16)
    # ax2.set_xlim([0,2])
    # ax2.set_ylim([0,120])
    ax2.legend(loc=0)
    plt.show()
