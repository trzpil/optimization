# plot a 2D graph for all damping without squeez film
# 2020-03-13 - First Version
# 2020-03-27 - Add Q_acoustic
# 2020-06-10 - update for paper

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
import numpy as np
import time
import sys
import os
# liste des repertoires ou se trouvent les fichiers
path = os.getcwd()
titre = sys.argv[0].replace('PV2I-',"").replace(".py","").replace(path+"/","")
print ("file name :",titre)

from fluid_material import fluid_material as fm
from geometry import cantilever as ct
from geometry import quartz_tuning_fork as qtf

#################################################
###   PARAMs                   
#################################################
start = time.time()

Air = fm.Fluid(name='air')
Si110 = fm.Material(name='silicon_100')

nbr_points = 200

Width = np.linspace(10e-6,10e-3,nbr_points) #[m]
Thickness = np.linspace(10e-6,10e-4,nbr_points) #[m]
Gap = 500e-6 #[m]
Freq = 15e3 #[Hz]

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
        'XiS':np.zeros((np.size(W,0),np.size(W,1)))}

#################################################
###   FUNCTIONs                   
#################################################

def Compute_cant(w,t):
    # w the width
    # t the thickness
    dico = {}
    Cant = ct.Cantilever(material=Si110,fluid=Air,dico={'name':'cant',
            'length':None,
            'width':w,
            'thickness':t,
            'gap':Gap,
            'freq':Freq,
            'Qvac':None})
    Cant.length_unkonw()
    dico['Qvis_sq'] = Cant.Q_viscous_sq()
    dico['Qsup'] = Cant.Q_support()
    dico['Qted'] = Cant.Q_thermoelastic()
    dico['Qaco'] = Cant.Q_acoustic()
    Cant.Qtot = (1/dico['Qvis_sq'] + 1/dico['Qsup'] + 1/dico['Qted'] + 1/dico['Qaco'])**(-1)
    dico['Qtot'] = Cant.Qtot
    # dico['Length'] = Cant.length
    # dico['Freq'] = Cant.f0()
    # dico['Surface'] = Cant.surface()
    # dico['Mass_eff'] = Cant.effective_mass()
    # dico['Xi'] = Cant.Xi0()
    # dico['XiS'] = dico['Xi']*dico['Surface']
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


cmap = 'Purples'

nbr_lvl = 25 # nrb of level for le colorbar

gridsize = (1, 1)
fig, _ = plt.subplots()

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
cs1 = ax1.contourf(T*1e6, W*1e6, data['Qtot'], cmap=cmap, levels = nbr_lvl)
ax1.contour(cs1, colors='k')
ax1.set_title(r'Q$_{total}$ (without squeeze) for f=%.1f kHz'%(Freq*1e-3))
ax1.set_ylabel('width (um)', labelpad=0)
ax1.set_xlabel('thickness (um)', labelpad=0)
ax1.set_xscale('log')
ax1.set_yscale('log')
cb = plt.colorbar(cs1, ax=ax1) # grab the Colorbar instance
for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)


# save figure
fig.savefig(titre+'_{:.0f}kHz.png'.format(Freq*1e-3))

plt.show()