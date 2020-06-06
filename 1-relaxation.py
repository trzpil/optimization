'''
Plot the acoustic pressure for CH4 as a function of the resonance frequency
2020-05-29 - first version
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

# liste des repertoires ou se trouvent les fichiers
path = os.getcwd()
titre = sys.argv[0].replace('PV2I-',"").replace(".py","").replace(path+"/","")
print ("file name :",titre)

from fluid_material import fluid_material as fm
from acoustic_pressure import acoustic_pressure as ap

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

AcPress = ap.AcousticPressure(fluid=Air,target_gas=Ch4,
        laser={'freq_mod':None,
            'laser_power':10e-3,
            'wavelength':1.65e-6,
            'waist':100e-6,
            'Rayleigh_length':None,
            'position':None})
AcPress.info()

zL = 150e-6 # [m] distance between laser beam and resonator
Frequecy = np.linspace(1e3,60e3,500) # frequency
Pressure_10000ppmv = [] # acoustic pressure 
for i in Frequecy:
    AcPress.freq_mod = i
    temp = AcPress.P(0,0,zL,0)
    Pressure_10000ppmv.append(temp)

Ch4.concentration = 0.005
Ch4.update()
Pressure_5000ppmv = [] # acoustic pressure 
for i in Frequecy:
    AcPress.freq_mod = i
    temp = AcPress.P(0,0,zL,0)
    Pressure_5000ppmv.append(temp)

######## MAX ##########
Pmax = max(np.abs(Pressure_10000ppmv))
index = np.where(np.abs(Pressure_10000ppmv) == Pmax)
print('\n')
print('Pressure is max for f = {} kHz'.format(Frequecy[index]*1e-3))

######## GRAPH ##########
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
gridsize = (1, 1)
fig, _ = plt.subplots()

number_color = 20
cmap = plt.get_cmap('tab20c')
colors = [cmap(i) for i in np.linspace(0, 1, number_color)]

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1.plot(Frequecy*1e-3, np.abs(Pressure_10000ppmv)*1e6, color=colors[0],
  linewidth=2,label=r'1% of CH$_4$ in N$_2$')
ax1.plot(Frequecy*1e-3, np.abs(Pressure_5000ppmv)*1e6, color=colors[2],
  linewidth=2,label=r'0.5% of CH$_4$ in N$_2$')
ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax1.tick_params(axis='both', which='major')
ax1.set_xlabel('modulation frequency (kHz)', labelpad=5)
ax1.set_ylabel(r'pressure amplitude ($\mu$Pa)', labelpad=5)
ax1.set_title(r'CH$_4$ - acoustic pressure')
ax1.legend(loc=0)

fig.savefig(titre+'.png')

plt.show()