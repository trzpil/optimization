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
titre = sys.argv[0].replace(".py","").replace(path+"/","")
print ("file name :",titre)

from fluid_material import fluid_material as fm
from geometry import cantilever as ct
from acoustic_pressure import acoustic_pressure as ap
from transduction import capacitance as cp 

#################################################
###   PARAMs                   
#################################################

nbr_points = 50

# Define resonator parameters
Width = np.linspace(1e-6,80e-6,nbr_points) #[m]
Thickness = np.linspace(1e-5,150e-6,nbr_points) #[m]
Gap = 4e-6 #[m]
Freq = 15e3 #[Hz]

#################################################
###   DEFINE CLASS INSTANCE                  
#################################################

# External fluid
Air = fm.Fluid(name='air')
Air.info()

# Targer gas
Ch4 = fm.GasSpecie(dico={
    'name':'CH4',
    'abs wavelength [m]':1.65e-6,
    'pressure [Pa]':Air.pressure,
    'temperature [K]':Air.temperature,
    'cross section [cm2/molecule]':1.55e-20,
    'relaxation time [s]':2.10e-5,
    'concentration [Nbr_mlc/Nbr_tot]':0.01})
Ch4.info()

# Cantilever material
Si100 = fm.Material(name='silicon_100')
Si100.info()

# Define cantilever
Cant = ct.Cantilever(material=Si100,fluid=Air,dico={'name':'cant',
        'length':None,
        'width':None,
        'thickness':None,
        'gap':Gap,
        'freq':Freq,
        'Qvac':None})
Cant.info()

# Acoustic pressure
AcPress = ap.AcousticPressure(fluid=Air,beam=Cant,target_gas=Ch4,
        laser={'freq_mod':Freq,
            'laser_power':10e-3,
            'wavelength':1.65e-6,
            'waist':100e-6,
            'Rayleigh_length':None,
            'position':[0,0,150e-6]}) # Define waist position of the laser
AcPress.info()

# Capacitance
Capa = cp.Capacitance(cantilever=Cant, acoustic=AcPress, Vdc=1)

#################################################
###   COMPUTE                    
#################################################
start = time.time()

W, T = np.meshgrid(Width, Thickness)

data = {'Thickness':np.zeros((np.size(W,0),np.size(W,1))),
        'Width':np.zeros((np.size(W,0),np.size(W,1))),
        'Qvis_sq':np.zeros((np.size(W,0),np.size(W,1))), # Q viscous with squeeez film 
        'Qvis_nosq':np.zeros((np.size(W,0),np.size(W,1))), # Q viscous without squeeze film
        'Qsup':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qted':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qaco':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qtot':np.zeros((np.size(W,0),np.size(W,1))), 
        'Qtot_nosq':np.zeros((np.size(W,0),np.size(W,1))),  # Qtotal without squeeze film
        'Length':np.zeros((np.size(W,0),np.size(W,1))),
        'Freq':np.zeros((np.size(W,0),np.size(W,1))), #resonnance frequency
        'Surface':np.zeros((np.size(W,0),np.size(W,1))), 
        'Mass_eff':np.zeros((np.size(W,0),np.size(W,1))), #effective mass
        'Xi':np.zeros((np.size(W,0),np.size(W,1))), # susceptibility
        'XiS':np.zeros((np.size(W,0),np.size(W,1))), 
        'W':np.zeros((np.size(W,0),np.size(W,1))), # displacement
        'Force':np.zeros((np.size(W,0),np.size(W,1))), # photo-acoustic force
        'C0':np.zeros((np.size(W,0),np.size(W,1))), # capacitance without any force
        'dVout':np.zeros((np.size(W,0),np.size(W,1)))} # amplitude of the output tension


def Compute_cant(w,t):
    # w the width
    # t the thickness
    dico = {}
    Cant.width = w
    Cant.thickness = t
    Cant.length_unkonw()

    AcPress.xL = 0.725*Cant.length # [m] value obtain from plot_beam_position.py

    # Fill the dictionary
    dico['Thickness'] = Cant.thickness
    dico['Width'] = Cant.width
    dico['Qvis_sq'] = Cant.Q_viscous_sq() 
    dico['Qvis_nosq'] = Cant.Q_viscous_nosq() # Q viscous without squeeze film
    dico['Qsup'] = Cant.Q_support()
    dico['Qted'] = Cant.Q_thermoelastic()
    dico['Qaco'] = Cant.Q_acoustic()
    Cant.Qtot = (1/dico['Qvis_sq'] + 1/dico['Qsup'] + 1/dico['Qted'] + 1/dico['Qaco'])**(-1)
    dico['Qtot'] = Cant.Qtot
    # dico['Qtot_nosq'] = (1/dico['Qvis_nosq'] + 1/dico['Qsup'] + 1/dico['Qted'] + 1/dico['Qaco'])**(-1)
    # dico['Length'] = Cant.length
    # dico['Freq'] = Cant.f0()
    # dico['Surface'] = Cant.surface()
    # dico['Mass_eff'] = Cant.effective_mass()
    # dico['Xi'] = Cant.Xi0()
    # dico['XiS'] = dico['Xi']*dico['Surface']
    # dico['Force'] = AcPress.Force()
    # dico['W'] = dico['Xi']*dico['Force']
    dico['C0'] = Capa.Czero()
    dico['dVout'] = Capa.dVout()
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


cmap = 'Greens'

nbr_lvl = 25 # nrb of level for le colorbar

gridsize = (1, 1)
fig, _ = plt.subplots()

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
cs1 = ax1.contourf(T*1e6, W*1e6, data['dVout']*1e6, cmap=cmap, levels = nbr_lvl)
ax1.contour(cs1, colors='k')
ax1.set_title(r'V$_{out}$ amplitude ($\mu$V) for f=%.0f kHz'%(Freq*1e-3))
ax1.set_ylabel('width (um)', labelpad=0)
ax1.set_xlabel('thickness (um)', labelpad=0)
# ax1.set_xscale('log')
# ax1.set_yscale('log')
cb = plt.colorbar(cs1, ax=ax1) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)


# save figure
fig.savefig(titre+'_{:.0f}kHz_{:.0f}um.png'.format(Freq*1e-3,Cant.gap*1e6))

plt.show()