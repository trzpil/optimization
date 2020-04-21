'''calculate the acoustic froce for different positions of the laser 
in regard to the cantilever
#### REF :
[1] Petra, N., Zweck, J., Kosterev, A. A., Minkoff, S. E., & Thomazy, D. (2009). 
Theoretical analysis of a quartz-enhanced photoacoustic spectroscopy sensor. 
Applied Physics B: Lasers and Optics, 94(4), 673–680. 
https://doi.org/10.1007/s00340-009-3379-1
/! see also pdf and matlab file from her student
#### VERSION
2020-04-17 - First version
2020-04-20 - update to take into acount the width of the laser
'''
import pandas #to write csv file
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
e = 40e-6 # thickness [m]
l =  100e-6 # width [m]
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


'''
Axis definition :
X is along the length of the cantilever
Y is along the width of the cantilever
Z is along the thickness of the cantilever
'''

# Define waist position of the laser
xL = 0 # [m]
yL = 0 # [m]
zL = 150e-6 # [m]

# Acoustic pressure
AcPress = ap.AcousticPressure(fluid=Air,
            dico={'freq_mod':Cant.freq,
                'laser_power':61.7e-3,
                'wavelength':1.53e-6,
                'waist':15e-6,
                'Rayleigh_length':None,
                'alpha_eff':1.31e-2})
AcPress.info()

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
    # ie along y-axis
    time = 0
    b = x-xL
    a = y-yL
    c = zL
    P1 = np.abs(AcPress.P(a,b,c,time))
    P2 = np.abs(AcPress.P(a,b,c+Cant.thickness,time))
    return P1-P2

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

def Force_density_perpendicular(x,y):
    return DeltaP_perpendicular(x,y)*phi(x)

def Force_density_parallel(x,y):
    return DeltaP_parallel(x,y)*phi(x)


#################################################
###   Compute force variation                  
#################################################

Force = {'perp':{'xL':[],'yL':[],'zL':[]}
        ,'para':{'xL':[],'yL':[],'zL':[]}}

nbr_points = 101

### xL variation 
tab_xL = np.linspace(0, 10*Cant.length, nbr_points)
yL = 0
for i in tab_xL:
    xL = i
    temp, _ = dblquad(Force_density_perpendicular, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['perp']['xL'].append(temp)
    temp, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['para']['xL'].append(temp)

# find where the perpendicular force is max
max_force_xL = max(Force['perp']['xL'])
index = Force['perp']['xL'].index(max_force_xL)
xL_perp_max = tab_xL[index] 
print('Force perpendicualr is max for xL=%.2fL'%(xL_perp_max/Cant.length))
# find where the parallel force is max
max_force_xL = max(Force['para']['xL'])
index = Force['para']['xL'].index(max_force_xL)
xL_para_max = tab_xL[index] 
print('Force parallel is max for xL=%.2fL'%(xL_para_max/Cant.length))

### yL variation 
tab_yL = np.linspace(-100*Cant.width, +100*Cant.width, nbr_points)
for i in tab_yL:
    yL = i
    xL = xL_perp_max
    temp, _ = dblquad(Force_density_perpendicular, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['perp']['yL'].append(temp)
    xL = xL_para_max
    temp, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['para']['yL'].append(temp)

# find where the perpendicular force is max
max_force_yL = max(Force['perp']['yL'])
index = Force['perp']['yL'].index(max_force_yL)
yL_perp_max = tab_yL[index] 
print('Force perpendicualr is max for yL=%.2fL'%(yL_perp_max/Cant.width))
# find where the parallel force is max
max_force_yL = max(Force['para']['yL'])
index = Force['para']['yL'].index(max_force_yL)
yL_para_max = tab_yL[index] 
print('Force parallel is max for yL=%.2fL'%(yL_para_max/Cant.width))


# zL variation 
tab_zL = np.linspace(zL, 2e-3, nbr_points)
for i in tab_zL:
    zL = i
    xL = xL_perp_max
    yL = yL_perp_max
    temp, _ = dblquad(Force_density_perpendicular, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['perp']['zL'].append(temp)
    xL = xL_para_max
    yL = yL_para_max
    temp, _ = dblquad(Force_density_parallel, -Cant.width/2, Cant.width/2, lambda toto: 0, lambda toto: Cant.length)
    Force['para']['zL'].append(temp)

#################################################
###   SAVE DATA                  
#################################################
info_cant = 'freq=%.1fkHz length=%.0fmm width=%.0fum thickness=%.0fum'%(Cant.freq*1e-3, Cant.length*1e3, Cant.width*1e6, Cant.thickness*1e6)
info_laser = 'mod=%.1fkHz lambda=%.2fum waist=%.0fum Rayleigh=%.1fmm alpha_eff=%fcm-1'%(AcPress.freq_mod*1e-3,AcPress.wavelength*1e6,AcPress.waist*1e6,AcPress.Rayleigh*1e3,AcPress.alpha_eff*1e-2)


donnees = {'x':tab_xL,
        'para_xL':Force['para']['xL'],
        'perp_xL':Force['perp']['xL'],
        'y':tab_yL,
        'para_yL':Force['para']['yL'],
        'perp_yL':Force['perp']['yL'],
        'z':tab_zL,
        'para_zL':Force['para']['zL'],
        'perp_zL':Force['perp']['zL'],
        }
data = pandas.DataFrame(donnees)
data.to_csv('beam_position - %s.csv'%(info_cant), sep=';', encoding='utf-8') 


#################################################
###   GRAPH                  
#################################################


gridsize = (2, 2)
fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
fig.suptitle('Beam Position (x-axis along cantilever length and Y-axis along the width) \n CANTILEVER : %s \n LASER : %s'%(info_cant,info_laser), fontsize=16)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1.plot(tab_xL/Cant.length, Force['para']['xL'],linewidth=2 ,
    label='Laser beam along x-axis (parallel) (max x/length=%.2f)'%(xL_para_max/Cant.length))
ax1.plot(tab_xL/Cant.length, Force['perp']['xL'],linewidth=2 
    ,label='Laser beam along y-axis (perpendicular) (max x/length=%.2f)'%(xL_perp_max/Cant.length))
ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_xlabel('Laser beam posotion on x-axis (x/length)', labelpad=5, fontsize=14)
ax1.set_ylabel('Acoustic force (N)', labelpad=5, fontsize=14)
# ax1.set_xlim([0,5])
ax1.set_ylim([0,2.7e-13])
ax1.legend(loc=0, prop={'size': 10})

ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax2.plot(tab_yL/Cant.width, Force['para']['yL'],linewidth=2 
    ,label='Laser beam along x-axis (parallel) (max y/width=%.0f)'%(yL_para_max/Cant.width))
ax2.plot(tab_yL/Cant.width, Force['perp']['yL'],linewidth=2 
    ,label='Laser beam along y-axis (perpendicular) (max y/width=%.0f)'%(yL_perp_max/Cant.width))
ax2.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.set_xlabel('Laser beam posotion on y-axis (x/width)', labelpad=5, fontsize=14)
ax2.set_ylabel('Acoustic force (N)', labelpad=5, fontsize=14)
# ax2.set_xlim([0,5])
ax2.set_ylim([0,2.7e-13])
ax2.legend(loc=0, prop={'size': 10})

ax3 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)
ax3.plot(tab_zL*1e3, Force['para']['zL'],linewidth=2 ,
    label='Laser beam along x-axis (parallel)')
ax3.plot(tab_zL*1e3, Force['perp']['zL'],linewidth=2 ,
    label='Laser beam along y-axis (perpendicular)')
ax3.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.set_xlabel('Laser beam posotion on z-axis (mm)', labelpad=5, fontsize=14)
ax3.set_ylabel('Acoustic force (N)', labelpad=5, fontsize=14)
ax2.set_ylim([0,2.7e-13])
ax3.legend(loc=0, prop={'size': 10})

fig.savefig('beam_position - %s .png'%(info_cant), dpi=fig.dpi*2)
plt.show()



