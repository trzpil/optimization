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
titre = sys.argv[0].replace('PV2I-',"").replace(".py","").replace(path+"/","")
print ("file name :",titre)

from fluid_material import fluid_material as fm
from geometry import cantilever as ct
from acoustic_pressure import acoustic_pressure as ap

#################################################
###   Calculation                  
#################################################

# Define fluid and material
Air = fm.Fluid(name='air')

Ch4 = fm.GasSpecie(dico={
    'name':'CH4',
    'abs wavelength [m]':1.65e-6,
    'pressure [Pa]':Air.pressure,
    'temperature [K]':Air.temperature,
    'cross section [cm2/molecule]':1.55e-20,
    'relaxation time [s]':2.10e-5,
    'concentration [Nbr_mlc/Nbr_tot]':0.01})
Ch4.info()

Si110 = fm.Material(name='silicon_110')

# Define cantilever
e = 100e-6 # thickness [m]
l =  500e-6 # width [m]
freq = 15e3 # resonance frequency [Hz]
Cant = ct.Cantilever(material=Si110,fluid=Air,dico={'name':'cantilever',
            'length':None,
            'width':l,
            'thickness':e,
            'gap':None,
            'freq':freq,
            'Qvac':None})
Cant.length_unkonw()
print('Cantilever length is : %f mm'%(Cant.length*1e3))
Cant.info()




'''
Axis definition :
X is along the length of the cantilever
Y is along the width of the cantilever
Z is along the thickness of the cantilever
'''

# Define waist position of the laser
xL = 0.725*Cant.length # [m] value obtain from plot_beam_position.py
yL = 0 # [m]
zL = 150e-6 # [m]

# Acoustic pressure with values from [1]
AcPress = ap.AcousticPressure(fluid=Air,target_gas=Ch4,
        laser={'freq_mod':Cant.freq,
            'laser_power':10e-3,
            'wavelength':1.65e-6,
            'waist':100e-6,
            'Rayleigh_length':None})
AcPress.info()

def phi(x):
    # return shape of the cantilever
    # with normalization with Phi_max =1
    # for the first mode phi is max for x=L 
    # where L is the length of the cantilever
    return Cant.phi(x)/Cant.phi(Cant.length)

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

def Force_density_parallel(x,y):
    return DeltaP_parallel(x,y)*phi(x)


#################################################
###   Compute force variation                  
#################################################

start = time.time()

nbr_points = 20
width_list = np.linspace(5e-6,800e-6,nbr_points)
thickness_list = np.linspace(5e-6,200e-6,nbr_points)
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
        force, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
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
# fig.suptitle('Acoustic force (pN)')

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
cs1 = ax1.contourf(thickness_grid*1e6, width_grid*1e6, data_plot*1e12, cmap=cmap, levels = nbr_lvl)
ax1.contour(cs1, colors='k')
ax1.set_title(r'acoustic force (pN) with f=%.1f kHz'%(AcPress.freq_mod*1e-3))
ax1.set_ylabel('width (um)', labelpad=0)
ax1.set_xlabel('thickness (um)', labelpad=0)
# ax1.grid(c='k', ls='-', alpha=0.3)
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ax1.tick_params(axis='both', which='major', labelsize=20)
# ax1.tick_params(axis='both', which='minor', labelsize=18)

#plt.colorbar(cs1, ax=ax1)
# increase color bat font size
cb = plt.colorbar(cs1, ax=ax1) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)


# save figure
fig.savefig(titre+'_{:.0f}kHz.png'.format(AcPress.freq_mod*1e-3))

plt.show()

