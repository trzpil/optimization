'''
Calculate the acoustic froce for different cantilever geometries
#### REF :
[1] Petra, N., Zweck, J., Kosterev, A. A., Minkoff, S. E., & Thomazy, D. (2009). 
Theoretical analysis of a quartz-enhanced photoacoustic spectroscopy sensor. 
Applied Physics B: Lasers and Optics, 94(4), 673–680. 
https://doi.org/10.1007/s00340-009-3379-1
/! see also pdf and matlab file from her student
#### VERSION
2020-05-20 - First version
2020-05-31 - Add relaxation time
'''
import time
import pandas #to write csv file
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.integrate import dblquad  # for double integration
import sys
import os

# liste des repertoires ou se trouvent les fichiers
path = os.getcwd()
titre = sys.argv[0].replace(".py","").replace(path+"/","")
print ("file name :",titre)

from fluid_material import fluid_material as fm
from geometry import cantilever as ct
from acoustic_pressure import acoustic_pressure as ap

#################################################
###   Parameters                 
#################################################
freq = 10e3 # resonance frequency [Hz]

nbr_points = 20
width_list = np.linspace(5e-6,800e-6,nbr_points)
thickness_list = np.linspace(5e-6,200e-6,nbr_points)

#################################################
###   Define instance of class                   
#################################################

# Define fluid  
Air = fm.Fluid(name='air')
Air.info()

# Define target gas
Ch4 = fm.GasSpecie(dico={
    'name':'CH4',
    'abs wavelength [m]':1.65e-6,
    'pressure [Pa]':Air.pressure,
    'temperature [K]':Air.temperature,
    'cross section [cm2/molecule]':1.55e-20,
    'relaxation time [s]':2.10e-5,
    'concentration [Nbr_mlc/Nbr_tot]':0.01})
Ch4.info()

# Define fluid and material 
Si100 = fm.Material(name='silicon_100')
Si100.info()

# Define cantilever
Cant = ct.Cantilever(material=Si100,fluid=Air,dico={'name':'cantilever',
            'length':None,
            'width':None,
            'thickness':None,
            'gap':None,
            'freq':freq,
            'Qvac':None})
Cant.info()

# Acoustic pressure with values from [1]
AcPress = ap.AcousticPressure(fluid=Air,target_gas=Ch4,beam=Cant,
        laser={'freq_mod':Cant.freq,
            'laser_power':10e-3,
            'wavelength':1.65e-6,
            'waist':100e-6,
            'Rayleigh_length':None,
            'position':None})
AcPress.info()


#################################################
###   Compute force variation                  
#################################################

start = time.time()

width_grid, thickness_grid = np.meshgrid(width_list, thickness_list)

data_plot = np.zeros((np.size(width_grid,0),np.size(width_grid,1)))
data_save = {'width [m]':[],'thickness [m]':[],'force [N]':[]}

for i in range(0,np.size(width_grid,0)):
    for j in range(0,np.size(width_grid,1)):
        width = width_grid[i][j]
        thickness = thickness_grid[i][j]
        Cant.width = width
        Cant.thickness = thickness
        Cant.length_unkonw()
        # Define waist position of the laser
        xL = 0.725*Cant.length # [m] value obtain from plot_beam_position.py
        yL = 0 # [m]
        zL = 150e-6 # [m]
        AcPress.xL, AcPress.yL, AcPress.zL = [xL,yL,zL]
        force = AcPress.Force()
        data_plot[i][j] = force
        data_save['width [m]'].append(width)
        data_save['thickness [m]'].append(thickness)
        data_save['force [N]'].append(force)

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


cmap = 'Blues'

nbr_lvl = 25 # nrb of level for le colorbar

gridsize = (1, 1)
fig, _ = plt.subplots()

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
cs1 = ax1.contourf(thickness_grid*1e6, width_grid*1e6, data_plot*1e12, cmap=cmap, levels = nbr_lvl)
ax1.contour(cs1, colors='k')
ax1.set_title(r'acoustic force (pN) for f=%.0f kHz'%(AcPress.freq_mod*1e-3))
ax1.set_ylabel('width (um)', labelpad=0)
ax1.set_xlabel('thickness (um)', labelpad=0)

cb = plt.colorbar(cs1, ax=ax1) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)


# save figure
fig.savefig(titre+'_{:.0f}kHz.png'.format(AcPress.freq_mod*1e-3))

plt.show()

