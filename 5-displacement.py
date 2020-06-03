# plot a 2D of W(x)
# 2020-03-13 - First Version
# 2020-03-27 - Add Q_acoustic
# 2020-06-03 - update for paper

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.integrate import dblquad  # for double integration
import time
import sys
import os
# liste des repertoires ou se trouvent les fichiers
path = os.getcwd()
titre = sys.argv[0].replace('PV2I-',"").replace(".py","").replace(path+"/","")
print ("file name :",titre)

from fluid_material import fluid_material as fm
from geometry import cantilever as ct
from acoustic_pressure import acoustic_pressure as ap 

#################################################
###   PARAMs                   
#################################################

nbr_points = 100

# Define resonator parameters
Width = np.linspace(1e-6,80e-6,nbr_points) #[m]
Thickness = np.linspace(1e-5,150e-6,nbr_points) #[m]
Gap = 4e-6 #[m]
Freq = 15e3 #[Hz]

#################################################
###   COMPUTE                  
#################################################


start = time.time()

Air = fm.Fluid(name='air')
Air.info()

Ch4 = fm.GasSpecie(dico={
    'name':'CH4',
    'abs wavelength [m]':1.65e-6,
    'pressure [Pa]':Air.pressure,
    'temperature [K]':Air.temperature,
    'cross section [cm2/molecule]':1.55e-20,
    'relaxation time [s]':2.10e-5,
    'concentration [Nbr_mlc/Nbr_tot]':0.01})
Ch4.info()

Si100 = fm.Material(name='silicon_100')
Si100.info()

# Acoustic pressure with values from [1]
AcPress = ap.AcousticPressure(fluid=Air,target_gas=Ch4,
        laser={'freq_mod':Freq,
            'laser_power':10e-3,
            'wavelength':1.65e-6,
            'waist':100e-6,
            'Rayleigh_length':None})
AcPress.info()


W, T = np.meshgrid(Width, Thickness)

data = {'Qvis_sq':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qsup':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qted':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qaco':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qtot':np.zeros((np.size(W,0),np.size(W,1))), 
        'Length':np.zeros((np.size(W,0),np.size(W,1))),
        'Freq':np.zeros((np.size(W,0),np.size(W,1))), 
        'Surface':np.zeros((np.size(W,0),np.size(W,1))), 
        'Mass_eff':np.zeros((np.size(W,0),np.size(W,1))), 
        'Xi':np.zeros((np.size(W,0),np.size(W,1))), 
        'XiS':np.zeros((np.size(W,0),np.size(W,1))),
        'W':np.zeros((np.size(W,0),np.size(W,1))),
        'Force':np.zeros((np.size(W,0),np.size(W,1)))}

#################################################
###   FUNCTIONs                   
#################################################

def Compute_cant(w,t):
    # w the width
    # t the thickness
    dico = {}
    Cant = ct.Cantilever(material=Si100,fluid=Air,dico={'name':'cant',
            'length':None,
            'width':w,
            'thickness':t,
            'gap':Gap,
            'freq':Freq,
            'Qvac':None})
    Cant.length_unkonw()

    # Define waist position of the laser
    xL = 0.725*Cant.length # [m] value obtain from plot_beam_position.py
    yL = 0 # [m]
    zL = 150e-6 # [m]

    # Cantilver shape mode use for the acoustic force
    def phi(x):
        # return shape of the cantilever
        # with normalization with Phi_max =1
        # for the first mode phi is max for x=L 
        # where L is the length of the cantilever
        return Cant.phi(x)/Cant.phi(Cant.length)

    # Pressure appplied on the cnatilever
    def DeltaP_parallel(x,y):
        # return the difference of pressure 
        # beteween the top and back of the mechanicla resoantor
        # when the laser beam is parallel to the length of the cantilever
        time = 0
        a = x-xL
        b = y-yL
        c = zL
        P1 = np.abs(AcPress.P(a,b,c,time))
        P2 = np.abs(AcPress.P(a,b,c+Cant.thickness,time))
        return P1-P2

    # Force on the cantilever
    def Force_density_parallel(x,y):
        return DeltaP_parallel(x,y)*phi(x)

    # Fill the dictionary
    dico['Qvis_sq'] = Cant.Q_viscous_sq()
    dico['Qsup'] = Cant.Q_support()
    dico['Qted'] = Cant.Q_thermoelastic()
    dico['Qaco'] = Cant.Q_acoustic()
    Cant.Qtot = (1/dico['Qvis_sq'] + 1/dico['Qsup'] + 1/dico['Qted'] + 1/dico['Qaco'])**(-1)
    dico['Qtot'] = Cant.Qtot
    # dico['Length'] = Cant.length
    # dico['Freq'] = Cant.f0()
    dico['Surface'] = Cant.surface()
    # dico['Mass_eff'] = Cant.effective_mass()
    dico['Xi'] = Cant.Xi0()
    # dico['XiS'] = dico['Xi']*dico['Surface']
    dico['Force'], _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    dico['W'] = dico['Xi']*dico['Force']
    return dico


for i in range(0,np.size(W,0)):
    for j in range(0,np.size(W,1)):
        results = Compute_cant(W[i][j],T[i][j])
        for k in results:
            data[k][i][j] = [results[k]][0]

end = time.time()
print('TIMER :')
info_time = 'exec_time=%.0fs =%.0fmin'%((end-start),(end-start)/60)
print(info_time) 


#################################################
###   PLOT DATA                               
#################################################

# MATPLOTLIB PARAMETERS
mpl.rcParams.update(mpl.rcParamsDefault)
# mpl.rcParamsDefault.keys --> display all parameters available
params = {'legend.fontsize': 16,
          'figure.figsize': (10, 8),
          'axes.labelsize': 18,
          'axes.titlesize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'figure.subplot.left': 0.15,
          'figure.subplot.right': 0.85,
          'figure.subplot.bottom': 0.15,
          'figure.subplot.top': 0.85,
          'xtick.direction': 'in',
          'ytick.direction': 'in',
          'xtick.major.size': 5,
          'xtick.major.width': 1.3,
          'xtick.minor.size': 3,
          'xtick.minor.width': 1,
          'xtick.major.pad': 8,
          'ytick.major.pad': 8,
          'lines.linewidth': 1,
          'axes.grid': 'True',
          'grid.alpha': 0.5,
          'grid.color': '111111',
          'grid.linestyle': '--',
          'savefig.dpi': 300, }

mpl.rcParams.update(params)


cmap = 'Oranges'

nbr_lvl = 25 # nrb of level for le colorbar

gridsize = (1, 1)
fig, _ = plt.subplots()

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
cs1 = ax1.contourf(T*1e6, W*1e6, data['W']*1e12, cmap=cmap, levels = nbr_lvl)
ax1.contour(cs1, colors='k')
ax1.set_title(r'W(x) displacement (pm) for f=%.1f kHz'%(Freq*1e-3))
ax1.set_ylabel('width (um)', labelpad=0)
ax1.set_xlabel('thickness (um)', labelpad=0)
# ax1.set_xscale('log')
# ax1.set_yscale('log')
cb = plt.colorbar(cs1, ax=ax1) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)


# save figure
fig.savefig(titre+'_{:.0f}kHz.png'.format(Freq*1e-3))

plt.show()