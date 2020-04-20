# calculate the acoustic froce for different positions of the laser 
# in regard to the cantilever
##### REF :
# [1] Petra, N., Zweck, J., Kosterev, A. A., Minkoff, S. E., & Thomazy, D. (2009). 
# Theoretical analysis of a quartz-enhanced photoacoustic spectroscopy sensor. 
# Applied Physics B: Lasers and Optics, 94(4), 673–680. 
# https://doi.org/10.1007/s00340-009-3379-1
# /!\ see also pdf and matlab file from her student
##### VERSION
# 2020-04-17 - First version

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import dblquad  # for double integration
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
from geometry import cantilever as ct
import acoustic_pressure as ap

#################################################
###   Calculation                  
#################################################

# Define fluid and material
Air = fm.Fluid(name='air')
Si110 = fm.Material(name='silicon_110')

# Define cantilever
e = 10e-6 # thickness [m]
l =  10e-6 # width [m]
freq = 10e3 # resonance frequency [Hz]
Cant = ct.Cantilever(material=Si110,fluid=Air,dico={'name':'cantilever',
            'length':None,
            'width':l,
            'thickness':e,
            'gap':None,
            'freq':freq,
            'Qvac':None})
Cant.length_unkonw()
print('Cantilever length is : %f mm'%(Cant.length*1e3))


# Define laser position
xL = 0.7*Cant.length # [m]
yL = 0 # [m]
zL = 150e-6 # [m]

# Acoustic pressure
AcPress = ap.AcousticPressure(Air)

def phi(x):
    # return shape of the cantilever
    # with normalization with Phi_max =1
    # for the first mode phi is max for x=L 
    # where L is the length of the cantilever
    return Cant.phi(x)/Cant.phi(Cant.length)

def DeltaP_perpendicular(x,y):
    # return the difference of pressure 
    # beteween the top and back of the mechanicla resoantor
    # when the laser beam is perpendicular to the length of the cantilever
    time = 0
    r1 = np.sqrt((x-xL)**2+(zL)**2)
    P1 = np.abs(AcPress.P(r1,time))
    r2 = np.sqrt((x-xL)**2+(e+zL)**2)
    P2 = np.abs(AcPress.P(r2,time))
    return P1-P2

def DeltaP_parallel(x,y):
    # return the difference of pressure 
    # beteween the top and back of the mechanicla resoantor
    # when the laser beam is parallel to the length of the cantilever
    time = 0
    r1 = np.sqrt((y-yL)**2+(zL)**2)
    P1 = np.abs(AcPress.P(r1,time))
    r2 = np.sqrt((y-yL)**2+(e+zL)**2)
    P2 = np.abs(AcPress.P(r2,time))
    return P1-P2

def Force_density_perpendicular(x,y):
    return DeltaP_perpendicular(x,y)*phi(x)

def Force_density_parallel(x,y):
    return DeltaP_parallel(x,y)*phi(x)


#################################################
###   Compute force variation                  
#################################################

Force = {'perp':{'xL':[],'yL':[],'zL':[]}
        ,'para':{'xL':[],'yL':[],'zL':[]}}

# xL variation 
tab_xL = np.linspace(0, 1.2*Cant.length, 50)
yL = 0
for i in tab_xL:
    xL = i
    temp, _ = dblquad(Force_density_perpendicular, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['perp']['xL'].append(temp)
    temp, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['para']['xL'].append(temp)

# find where the perpendicular for is max
max_force_xL = max(Force['perp']['xL'])
index = Force['perp']['xL'].index(max_force_xL)
xL_max = tab_xL[index] 
print('Force perpendicualr is max for xL=%fL'%(xL_max/Cant.length))

# yL variation 
tab_yL = np.linspace(-10*Cant.width, +10*Cant.width, 50)
xL = xL_max
for i in tab_yL:
    yL = i
    temp, _ = dblquad(Force_density_perpendicular, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['perp']['yL'].append(temp)
    temp, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['para']['yL'].append(temp)

# zL variation 
tab_zL = np.linspace(150e-6, 2e-3, 50)
xL = xL_max
yL=0
for i in tab_zL:
    zL = i
    temp, _ = dblquad(Force_density_perpendicular, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['perp']['zL'].append(temp)
    temp, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['para']['zL'].append(temp)

#################################################
###   GRAPH                  
#################################################
info = 'freq=%.1fkHz length=%.0fmm width=%.0fum thickness=%.0fum'%(Cant.freq*1e-3, Cant.length*1e3, Cant.width*1e6, Cant.thickness*1e6)

gridsize = (1, 3)
fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
fig.suptitle('Beam Position (x-axis belong cantilever length) \n CANTILEVER : '+info, fontsize=20)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1.plot(tab_xL*1e3, Force['perp']['xL'],linewidth=2 ,label='Laser beam along y-axis (perpendicular to the cantilever)')
ax1.plot(tab_xL*1e3, Force['para']['xL'],linewidth=2 ,label='Laser beam along x-axis (parallel to the cantilever)')
ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_xlabel('xL - Laser beam posotion on x-axis (mm)', labelpad=5, fontsize=14)
ax1.set_ylabel('Acoustic force (N)', labelpad=5, fontsize=14)
# ax1.set_xlim([0,5])
ax1.set_ylim([0,1.3e-14])
ax1.legend(loc=0)

ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax2.plot(tab_yL*1e6, Force['perp']['yL'],linewidth=2 ,label='Laser beam along y-axis (perpendicular to the cantilever)')
ax2.plot(tab_yL*1e6, Force['para']['yL'],linewidth=2 ,label='Laser beam along x-axis (parallel to the cantilever)')
ax2.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.set_xlabel('yL - Laser beam posotion on y-axis (um)', labelpad=5, fontsize=14)
ax2.set_ylabel('Acoustic force (N)', labelpad=5, fontsize=14)
# ax2.set_xlim([0,5])
ax2.set_ylim([0,1.3e-14])
ax2.legend(loc=0)

ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)
ax3.plot(tab_zL*1e3, Force['perp']['zL'],linewidth=2 ,label='Laser beam along y-axis (perpendicular to the cantilever)')
ax3.plot(tab_zL*1e3, Force['para']['zL'],linewidth=2 ,label='Laser beam along x-axis (parallel to the cantilever)')
ax3.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.set_xlabel('zL - Laser beam posotion on z-axis (mm)', labelpad=5, fontsize=14)
ax3.set_ylabel('Acoustic force (N)', labelpad=5, fontsize=14)
ax3.set_ylim([0,1.3e-14])
ax3.legend(loc=0)

fig.savefig('beam_position.png', dpi=fig.dpi*2)
plt.show()

