# Capacitance and its variation from X the susceptibility
# [1] Mu, L., Zhang, W., He, C., Zhang, R., Song, J., & Xue, C. (2015). 
# Design and test of a capacitance detection circuit based on a transimpedance 
# amplifier. 
# Journal of Semiconductors, 36(7), 1–8. 
# https://doi.org/10.1088/1674-4926/36/7/075007
# [2] red note book page 116

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
from acoustic_pressure import acoustic_pressure as ap 

from scipy import constants 
eps_0 = constants.physical_constants["electric constant"][0] # [F m^-1] the electric constant (vacuum permittivity)
print('vacuum permittivity : ',eps_0)

class Capacitance():
    # define a class wich can return the capacitance and its variations
    def __init__(self,cantilever,acoustic,Vdc):
        self.cant = cantilever
        self.acoustic = acoustic
        self.Vdc = Vdc #tension of polarisation [V]

    def Czero(self):
        # return the value of the capacitor for parallel electrode of surface S
        # separate by a distance d
        # S : surface [m2]
        # d : dictance between electrode [m]
        d = self.cant.gap
        S = self.cant.surface()
        eps_r = self.cant.fluid.relative_permittivity
        C = eps_0*eps_r*S/d
        return C # [F]

    def Cpa(self):
        '''
        compute capacitance when a acoustic force is applied
        '''
        # capacitance without any force [F]
        C0 = self.Czero()  
        # appplies force [N]
        F = self.acoustic.Force()
        # susceptibility [m/N]
        Chi = self.cant.Xi0()
        # cantilver length and gap [m]
        L = self.cant.length
        d = self.cant.gap
        # cantilever displacement in the z direction [m]
        W = F*Chi
        # integration from 0 to Length of phi(x) the first cantilver shape mode 
        int_phi = 0.783/2*self.cant.length # cf Aoust's thesis page 299
        # capacitance
        Cmax = C0 + C0/d/L*W*int_phi
        Cmin = C0 - C0/d/L*W*int_phi
        return [Cmin,Cmax]

    def dVout(self):
        # ampltide of the output tension
        # capacitance without any force [F]
        C0 = self.Czero() 
        C_min, C_max = self.Cpa()
        Vdc = self.Vdc
        Vout = (C0/C_min-C0/C_max)*Vdc
        return Vout      

if __name__ == "__main__":

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
    Si = fm.Material(name='silicon_100')
    Si.info()

    # Cantilever
    Cant = ct.Cantilever(material=Si,fluid=Air,dico={'name':'cant',
                'length':1e-3,
                'width':150e-6,
                'thickness':10e-6,
                'gap':5e-6,
                'freq':None,
                'Qvac':100})
    Cant.Qtot = Cant.Qvac
    Cant.info()

    # Acoustic pressure
    AcPress = ap.AcousticPressure(fluid=Air,beam=Cant,target_gas=Ch4,
            laser={'freq_mod':Cant.freq,
                'laser_power':10e-3,
                'wavelength':1.65e-6,
                'waist':100e-6,
                'Rayleigh_length':None,
                'position':[0,0,150e-6]}) # Define waist position of the laser
    AcPress.info()

    # Capacitor
    Capa = Capacitance(cantilever=Cant, acoustic=AcPress, Vdc=1)
    print('C0    = {} fF'.format(Capa.Czero()*1e15))
    print('dVout = {} uV'.format(Capa.dVout()*1e6))


