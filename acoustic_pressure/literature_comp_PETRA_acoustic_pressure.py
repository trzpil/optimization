# acoustic pressure modelisation
# The goal of this code is to find which equation of [1] we should use
##### REF :
# [1] Petra, N., Zweck, J., Kosterev, A. A., Minkoff, S. E., & Thomazy, D. (2009). 
# Theoretical analysis of a quartz-enhanced photoacoustic spectroscopy sensor. 
# Applied Physics B: Lasers and Optics, 94(4), 673–680. 
# https://doi.org/10.1007/s00340-009-3379-1
# /!\ see also pdf and matlab file from her student
##### VERSION
# 2020-04-09 - first version
# 2020-04-17 - change name of the file 

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

Air = fm.Fluid(name='air')
print('Air heat capacity ratio : ',Air.Gamma)

# frequency
freq = 32761 #[Hz]
omega = 2*np.pi*freq #[1/s]
# laser beam width
sigma = 0.05e-3 #[m]
# laser power
Pl = 61.7e-3 #[W]
# effecrive losses
k_eff = 1.31e-2 #[m-1]
# air heat capacity ratio
Gamma = Air.Gamma
# amplitude of the acoustic source
W = -1*(Gamma-1)*k_eff/2*Pl/(2*np.pi*sigma**2)*omega
print('Amplitude of the acoustic source : %g N/m2/s2'%W)
# speed of sound
c = Air.speed_of_sound
print('Speed of sound : %f m/s'%c)

def gaussian(r):
    # laser beam profile
    return np.exp(-1*r**2/2./sigma**2)

def f1(r):
    # equation (8) from [1]
    k = omega/c
    A = np.pi*W/2./c**2
    def f1A(s):
        return s*y0(k*s)*gaussian(s)
    # First integration
    int_f1A, _ = quad(f1A, r, np.inf) 
    def f1B(s):
        return s*j0(k*s)*gaussian(s)
    # Second integration 
    int_f1B, _ = quad(f1B, 0., r)
    return -A*(int_f1A*j0(k*r)+int_f1B*y0(k*r))

def f2(r):
    # equation (9) from [1]
    # using the change of variables u = ks otherwise f2(r)=0 
    # (Why ? I don't hnow ???)
    k = omega/c
    A = np.pi*W/2./c**2/k**2
    def f2A(s):
        #return s*j0(k*s)*gaussian(s)
        return s*j0(s)*np.exp(-s**2/2./k**2/sigma**2) # s is u
    int_f2A, _ = quad(f2A, 0., np.inf)
    return A*int_f2A*j0(k*r)

def f1_minus_i_f2(r):
    # equation (13) from [1]    
    k = omega/c
    A = -np.pi*W/2./c**2/k**2
    def h(u):
        return u*j0(u)*np.exp(-u**2/2./k**2/sigma**2)
    int_h, _ = quad(h, 0., np.inf)
    return A*int_h*(y0(k*r)+1j*j0(k*r))


# Acouistic pressure
def P(r,t):
    # equation (7) with (8) and (9) from [1]
    pr = f1(r)-1j*f2(r)
    return pr*np.exp(1j*omega*t)

def P_bis(r,t):
    # equation (7) with (13) from [1]
    pr = f1_minus_i_f2(r)
    return pr*np.exp(1j*omega*t)

def P_simplified(r,t):
    # equation (14) from [1]
    k = omega/c
    A = (Gamma-1)*omega*k_eff*Pl/8/c**2
    part1 = j0(k*r)*np.cos(omega*t)
    part2 = y0(k*r)*np.sin(omega*t)
    return A*(part1+part2)

def amp_P_simple(r):
    # return the amplitude of P_simplified
    T = 1/freq
    time = np.linspace(0,3*T,100)
    wave = P_simplified(r,time)
    amp = 1./2.*(max(wave)-min(wave))
    return amp


#################################################
###   PLOT DATA                               
#################################################
# data
r = np.linspace(1e-9,2e-3,1000)
P = np.asarray([P(i,0) for i in r])
P_bis = np.asarray([P_bis(i,0) for i in r])
P_simplified = np.asarray([amp_P_simple(i) for i in r])

# color list
cmap = plt.get_cmap('tab10')
number = 10
colors = [cmap(i) for i in np.linspace(0, 1, number)]

# graph
gridsize = (1, 2)
fig, _ = plt.subplots(figsize=(20, 10), dpi=90)
fig.suptitle('PETRA Acoustic pressure', fontsize=20)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1.plot(r*1e3, np.abs(P)*1e3,linewidth=2, color=colors[0],label='equation (7)')
ax1.plot(r*1e3, np.abs(P_bis)*1e3,linewidth=2, color=colors[1],label='equation (7) with (13)')
ax1.plot(r*1e3, P_simplified*1e3,linewidth=2, color=colors[2],label='equation simplified (14)')
ax1.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_xlabel('Radial distance from laser beam (mm)', labelpad=5, fontsize=14)
ax1.set_ylabel('Pressure amplitude (mPa)', labelpad=5, fontsize=14)
# ax1.set_title('title',fontsize=16)
ax1.set_xlim([0,2])
ax1.set_ylim([0,0.2])
ax1.legend(loc=3)

ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax2.plot(r*1e3, np.angle(P,deg=True),linewidth=2, color=colors[0],label='equation (7)')
ax2.plot(r*1e3, np.angle(P_bis,deg=True),linewidth=2, color=colors[1],label='equation (7) with (13)')
ax2.grid(b=True, which='major', color='#D3D3D3', linestyle='-')
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.set_xlabel('Radial distance from laser beam (mm)', labelpad=5, fontsize=14)
ax2.set_ylabel('Pressure phase (degrees)', labelpad=5, fontsize=14)
# ax2.set_title('title',fontsize=16)
# ax2.set_xlim([0,2])
# ax2.set_ylim([0,120])
ax2.legend(loc=0)

fig.savefig('literature_comp_PETRA_acoustic_pressure.png', dpi=fig.dpi*2)
plt.show()
