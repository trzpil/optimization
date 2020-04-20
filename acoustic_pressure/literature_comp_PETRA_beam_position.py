# compute acoustic force whe the laser beam is perpendicular to the length
# of the QTF
# the goal is to take into acount the shape of the QTF (the phi function)
# and reproduce figure 3 of [1] 
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

Air = fm.Fluid(name='air')
Air.info()
Quartz = fm.Material(name='quartz')
QTF = ct.Cantilever(material=Quartz,fluid=Air,dico={'name':'QTF_Petra',
            'length':3.8e-3,
            'width':340e-6,
            'thickness':600e-6,
            'gap':300e-6,
            'freq':32.8e3,
            'Qvac':None})

AcPress = ap.AcousticPressure(Air)

# Laser position
xL = 0.7*QTF.length
zL = QTF.gap/2

# resonnator
e = QTF.thickness
l = QTF.width
L = QTF.length

def phi(x):
    # return shape of the cantilever
    # with normalization with Phi_max =1
    # for the first mode phi is max for x=L 
    return QTF.phi(x)/QTF.phi(L)

def DeltaP(x,y):
    # return the difference of pressure 
    # beteween the top and back of the mechanicla resoantor
    time = 0
    r1 = np.sqrt((x-xL)**2+(zL)**2)
    P1 = np.abs(AcPress.P(r1,time))
    r2 = np.sqrt((x-xL)**2+(e+zL)**2)
    P2 = np.abs(AcPress.P(r2,time))
    return P1-P2

def Force_density(x,y):
    return DeltaP(x,y)*phi(x)


Force = []
Position = np.linspace(0,5e-3,100)
for i in Position:
    xL = i
    temp, _ = dblquad(Force_density, -l/2, l/2, lambda toto: 0, lambda toto: L)
    Force.append(temp)

# normalize the Force
Force = np.array(Force)/max(Force)

#################################################
###   Get roman data                 
#################################################
Data = np.loadtxt(open('roman produit de convolution/data-ROMAN.csv', "rb"), 
    delimiter=';', skiprows=1, unpack=True) 
#print(Data)
#################################################
###   GRAPH                  
#################################################


gridsize = (1, 1)
fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
fig.suptitle('PETRA Beam Position', fontsize=20)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1.plot(Position*1e3, Force,linewidth=2 ,label='with phi normalized')
ax1.plot(Data[1], Data[2],linewidth=2 ,label='roman simulation')
ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_xlabel('Beam posotion (mm)', labelpad=5, fontsize=14)
ax1.set_ylabel('Normalized amplitude', labelpad=5, fontsize=14)
ax1.set_xlim([0,5])
ax1.set_ylim([0,1])
ax1.legend(loc=0)

fig.savefig('literature_comp_PETRA_beam_position.png', dpi=fig.dpi*2)
plt.show()

