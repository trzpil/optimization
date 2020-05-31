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
    def __init__(self,fluid,target_gas,laser):
        self.gas = target_gas 
        self.fluid = fluid
        # modualtion frequency
        self.freq_mod = laser['freq_mod'] #[Hz]        
        # laser power
        self.Pl = laser['laser_power'] #[W]
        # laser wavelength
        self.wavelength = laser['wavelength'] # [m]
        # laser waist (minimum value possible)
        # if we compare with [1] sigma = waist/2
        if laser['waist'] is None:
            self.waist = 2*self.wavelength/np.pi
        else:
            self.waist = laser['waist']
        # laser Rayleigh length
        if laser['Rayleigh_length'] is None:
            self.Rayleigh = np.pi*self.waist**2/self.wavelength
        else:
            self.Rayleigh = laser['Rayleigh_length'] 
        # effecrive losses
        self.alpha_eff = self.gas.absorption_coeff # [m-1]

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
        # effecrive losses the 0.7 come from the result obtain for Keff in [1]
        # /! k_eff is divide by 2 in W
        k_eff = self.alpha_eff*0.7
        # relaxation time [s]
        tau = self.gas.relaxation_time
        # amplitude of the acoustic source
        W = -1*(Gamma-1)*k_eff/2*Pl*1./np.sqrt(1+(omega*tau)**2)*2./(np.pi*wL**2)*omega
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

    Ch4 = fm.GasSpecie(dico={
        'name':'CH4',
        'abs wavelength [m]':1.65e-6,
        'pressure [Pa]':Air.pressure,
        'temperature [K]':Air.temperature,
        'cross section [cm2/molecule]':1.55e-20,
        'relaxation time [s]':15e-6,
        'concentration [Nbr_mlc/Nbr_tot]':0.01})
    ################################################
    # change absorption coeff and relaxation time  #
    # to compare with Petra results                #
    ################################################
    Ch4.absorption_coeff = 1.31e-2/0.7
    Ch4.relaxation_time = 1e-9
    Ch4.info()

    AcPress = AcousticPressure(fluid=Air,target_gas=Ch4,
            laser={'freq_mod':32.8e3,
                'laser_power':61.7e-3,
                'wavelength':1.53e-6,
                'waist':100e-6,
                'Rayleigh_length':None})
    AcPress.info()

    Pressure = []
    tab_z = np.linspace(1e-9,2e-3,100)
    for i in tab_z:
        temp = AcPress.P(0,0,i,0)
        Pressure.append(temp)       
    Pressure = np.array(Pressure)

    # graph
    gridsize = (1, 2)
    fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
    fig.suptitle('PETRA Acoustic pressure', fontsize=20)

    '''
    PETRA results for:
    laser={'freq_mod':32.8e3,
                'laser_power':61.7e-3,
                'wavelength':1.53e-6,
                'waist':100e-6,
                'Rayleigh_length':None,
                'alpha_eff':1.31e-2}
    '''

    PETRA_Z = np.array([1.00000000e-09, 2.02030101e-05, 4.04050202e-05, 6.06070303e-05
        , 8.08090404e-05, 1.01011051e-04, 1.21213061e-04, 1.41415071e-04
        , 1.61617081e-04, 1.81819091e-04, 2.02021101e-04, 2.22223111e-04
        , 2.42425121e-04, 2.62627131e-04, 2.82829141e-04, 3.03031152e-04
        , 3.23233162e-04, 3.43435172e-04, 3.63637182e-04, 3.83839192e-04
        , 4.04041202e-04, 4.24243212e-04, 4.44445222e-04, 4.64647232e-04
        , 4.84849242e-04, 5.05051253e-04, 5.25253263e-04, 5.45455273e-04
        , 5.65657283e-04, 5.85859293e-04, 6.06061303e-04, 6.26263313e-04
        , 6.46465323e-04, 6.66667333e-04, 6.86869343e-04, 7.07071354e-04
        , 7.27273364e-04, 7.47475374e-04, 7.67677384e-04, 7.87879394e-04
        , 8.08081404e-04, 8.28283414e-04, 8.48485424e-04, 8.68687434e-04
        , 8.88889444e-04, 9.09091455e-04, 9.29293465e-04, 9.49495475e-04
        , 9.69697485e-04, 9.89899495e-04, 1.01010151e-03, 1.03030352e-03
        , 1.05050553e-03, 1.07070754e-03, 1.09090955e-03, 1.11111156e-03
        , 1.13131357e-03, 1.15151558e-03, 1.17171759e-03, 1.19191960e-03
        , 1.21212161e-03, 1.23232362e-03, 1.25252563e-03, 1.27272764e-03
        , 1.29292965e-03, 1.31313166e-03, 1.33333367e-03, 1.35353568e-03
        , 1.37373769e-03, 1.39393970e-03, 1.41414171e-03, 1.43434372e-03
        , 1.45454573e-03, 1.47474774e-03, 1.49494975e-03, 1.51515176e-03
        , 1.53535377e-03, 1.55555578e-03, 1.57575779e-03, 1.59595980e-03
        , 1.61616181e-03, 1.63636382e-03, 1.65656583e-03, 1.67676784e-03
        , 1.69696985e-03, 1.71717186e-03, 1.73737387e-03, 1.75757588e-03
        , 1.77777789e-03, 1.79797990e-03, 1.81818191e-03, 1.83838392e-03
        , 1.85858593e-03, 1.87878794e-03, 1.89898995e-03, 1.91919196e-03
        , 1.93939397e-03, 1.95959598e-03, 1.97979799e-03, 2.00000000e-03])

    PETRA_Pressure = np.array([6.38457677e-04, 2.10949354e-04, 1.82432625e-04, 1.66050225e-04
        , 1.54601602e-04, 1.45841286e-04, 1.38772924e-04, 1.32866792e-04
        , 1.27807635e-04, 1.23392668e-04, 1.19483791e-04, 1.15982745e-04
        , 1.12817152e-04, 1.09932175e-04, 1.07285282e-04, 1.04842831e-04
        , 1.02577759e-04, 1.00467977e-04, 9.84952267e-05, 9.66442436e-05
        , 9.49021408e-05, 9.32579412e-05, 9.17022193e-05, 9.02268237e-05
        , 8.88246584e-05, 8.74895090e-05, 8.62159036e-05, 8.49989995e-05
        , 8.38344910e-05, 8.27185334e-05, 8.16476799e-05, 8.06188288e-05
        , 7.96291796e-05, 7.86761952e-05, 7.77575702e-05, 7.68712041e-05
        , 7.60151774e-05, 7.51877320e-05, 7.43872537e-05, 7.36122567e-05
        , 7.28613712e-05, 7.21333309e-05, 7.14269636e-05, 7.07411821e-05
        , 7.00749759e-05, 6.94274045e-05, 6.87975914e-05, 6.81847180e-05
        , 6.75880189e-05, 6.70067777e-05, 6.64403224e-05, 6.58880225e-05
        , 6.53492852e-05, 6.48235526e-05, 6.43102994e-05, 6.38090299e-05
        , 6.33192763e-05, 6.28405963e-05, 6.23725718e-05, 6.19148067e-05
        , 6.14669257e-05, 6.10285727e-05, 6.05994098e-05, 6.01791159e-05
        , 5.97673855e-05, 5.93639281e-05, 5.89684670e-05, 5.85807384e-05
        , 5.82004908e-05, 5.78274842e-05, 5.74614894e-05, 5.71022873e-05
        , 5.67496686e-05, 5.64034328e-05, 5.60633880e-05, 5.57293507e-05
        , 5.54011446e-05, 5.50786010e-05, 5.47615578e-05, 5.44498595e-05
        , 5.41433568e-05, 5.38419062e-05, 5.35453697e-05, 5.32536146e-05
        , 5.29665133e-05, 5.26839429e-05, 5.24057851e-05, 5.21319257e-05
        , 5.18622550e-05, 5.15966670e-05, 5.13350594e-05, 5.10773335e-05
        , 5.08233943e-05, 5.05731496e-05, 5.03265107e-05, 5.00833918e-05
        , 4.98437097e-05, 4.96073844e-05, 4.93743381e-05, 4.91444956e-05])

    '''
    PLOT
    '''

    ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
    ax1.plot(tab_z*1e3, np.abs(Pressure)*1e3,linewidth=2,label='equation (7) with (13)')
    ax1.plot(PETRA_Z*1e3, PETRA_Pressure*1e3,linewidth=1,label='Petra results', color='k', dashes=[12, 10])
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
