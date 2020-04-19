# compute hydrodynamic function to estimate viscous damping
# REF [1]
# Aoust, G., Levy, R., Bourgeteau, B., & Le Traon, O. (2015). 
# Viscous damping on flexural mechanical resonators. 
# Sensors and Actuators, A: Physical, 230, 126–135. 
# https://doi.org/10.1016/j.sna.2015.04.004
# REF [2]
# Sader, J. E. (1998). 
# Frequency response of cantilever beams immersed in viscous fluids with 
# applications to the atomic force microscope. 
# Journal of Applied Physics, 84(1), 64–76. 
# https://doi.org/10.1063/1.368002
# REF [3]
# Green, C., & Sader, J. (2005).
# Frequency response of cantilever beams immersed in viscous fluids.
# Journal of Applied Physics.
# Retrieved from http://link.aip.org/link/?jap/98/114913
# REF [4]
# Cox, R., Josse, F., Heinrich, S. M., Brand, O., & Dufour, I. (2012).
# Characteristics of laterally vibrating resonant microcantilevers in viscous 
# liquid media. 
# Journal of Applied Physics, 111(1). 
# https://doi.org/10.1063/1.3674278


# 2020-02-05 add GammaVacuum
# 2020-02-18 modified calculation of density and dyna visocity in fluid class
# 2020-02-19 add quartz / correction squeeze
# 2020-02-20 add mass in Beam class
# 2020-03-04 add class cantilever
# 2020-03-06 - Import fluid and material from an other directory 
# 2020-03-15 - add Hydro function from ref[3]  take into account edge and thickness effect 
# 2020-03-20 - add auto calcul of Reynolds number when is needed

import numpy as np
from scipy.special import kv # modified Bessel functions of the second kind

import sys
import os
# get current working directory
path = os.path.abspath(os.getcwd())
index = path[::-1].find('/')
path = path[:len(path)-index-1]
print('current working directory : ',path)
# set the path of the module I want
sys.path.append(path) 

from fluid_material import fluid_material as fm
from geometry import beam as be



class HydrodynamicFunction():
    # The hydrodynamic function, is the total hydrody-
    # namic force per unit length normalized to the amount of
    # force per unit length it would take to excite fluid occupy-
    # ing a circular cylindrical volume with a diameter equal to
    # the microcantilever’s width to the same velocity as the microcantilever [4]
    def __init__(self,fluid,beam):
        self.fluid = fluid # Fluid object
        self.beam = beam # Beam object
        self.Re = self.ReynoldsNumber()
        self.Gfb = self.GammaFrontBack()
        self.Gup = self.GammaUpDown()
        if self.beam.Qvac is not None:
            self.Gvac = self.GammaVacuum()
        self.Sader = self.SaderCoefficient()
        if self.beam.gap is not None:
            self.Gsqueeze = self.GammaSqueeze()
        # use the superposition principle to retrieved 
        # the expression of the hydrodynamic function
        #self.Ghydro = self.Sader*self.Gfb + self.Gup
        # quality factor extract from Ghydro
        #self.Q = self.QualityFactor()

    def GammaFrontBack(self):
        # It is the the exact analytical result for a beam with a circular 
        # cross section
        # An analytical solution exists in terms of Meijer G-functions 
        # but is rather cumbersome and has no direct expression [3]
        Re = self.ReynoldsNumber()
        X = -1j*np.sqrt(1j*Re)
        # K are modified Bessel functions of the second kind
        K0 = kv(0,X)
        K1 = kv(1,X)
        return 1+4*1j*K1/K0/np.sqrt(1j*Re)

    def GammaUpDown(self):
        # See ref [1]
        # according to ref[4] this expression doesn't take into account 
        # edge and thickness effect
        Re = self.Re
        e = self.beam.thickness
        l = self.beam.width
        num = 2*np.sqrt(2)*e
        den = np.pi*l*np.sqrt(Re)
        return num/den*(1+1j)

    def GammaUpDown_COX(self,num_correction=True):
        # From ref[4] to take into account 
        # edge and thickness effect
        Re = self.Re
        h = self.beam.thickness
        b = self.beam.width
        if num_correction is True: 
            Cr = 1.658*(h/b)**1.83*np.sqrt(Re)+3.08*(h/b)**0.85+1
            Ci = (2.56-1.321*(h/b))*1/np.sqrt(Re)+3.108*(h/b)**0.85+1
        else:
            Cr=1
            Ci=1
        return 4/(np.pi*np.sqrt(2*Re))*(Cr+1j*Ci)

    def GammaSqueeze(self):
        omega = 2*np.pi*self.beam.freq # the radial frequency
        l = self.beam.width
        d = self.beam.gap # gap of air for the squeeze film
        mu = self.fluid.dyna_viscosity
        rho = self.fluid.density
        P = self.fluid.pressure
        # Squeeze number
        sigma = 12*mu*omega*l**2/P/d**2
        # Langlois's function
        X = np.sqrt(sigma/2)
        fe = 1 - 1/X * (np.sinh(X)+np.sin(X))/(np.cosh(X)+np.cos(X))
        fd = 1/X * (np.sinh(X)-np.sin(X))/(np.cosh(X)+np.cos(X))
        # Gamma squeeze film
        G = -4*P/np.pi/d/l/rho/omega**2*(fe-1j*fd)
        return G

    def GammaVacuum(self):
        # GammaVacuum is taken into account in order to maintain correct 
        # asymptotic behaviour when the viscous damping factor disappears 
        # (e.g. at very low pressure). An effective contribution is thus added
        # to GammaHydro the total hydrodynamic function.
        # [Aoust's thesis page 89]
        omega = 2*np.pi*self.beam.freq # the radial frequency
        l = self.beam.width
        e = self.beam.thickness
        rho_b = self.beam.density
        rho_f = self.fluid.density
        Q = self.beam.Qvac
        X = 4*rho_b*e/np.pi/rho_f/l
        return 1j*X/Q

    def QualityFactor(self):
        omega = 2*np.pi*self.beam.freq # the radial frequency
        l = self.beam.width
        e = self.beam.thickness
        rho_b = self.beam.density
        rho_f = self.fluid.density
        G = self.Ghydro # the total hydrodynamic function
        Gr = np.real(G) # real part of the total hydrodynamic function
        Gi = np.imag(G) # imaginary part of the total hydrodynamic function
        X = 4*rho_b*e/np.pi/rho_f/l
        # The expression of Q can be simplified as for most of the case 
        # X >> GammaFB
        # print('X : ',X)
        Q = (X+Gr)/Gi
        return Q

    def QualityFactor_COX(self):
        # From ref[4] to avoid approximation on Q
        omega = 2*np.pi*self.beam.freq # the radial frequency
        l = self.beam.width
        e = self.beam.thickness
        rho_b = self.beam.density
        rho_f = self.fluid.density
        G = self.Ghydro # the total hydrodynamic function
        GammaR = np.real(G) # real part of the total hydrodynamic function
        GammaI = np.imag(G) # imaginary part of the total hydrodynamic function
        X = 4/np.pi*rho_f*l**2
        g1 = X*GammaI*omega
        g2 = X*GammaR
        Q = (2*(1-np.sqrt(1-g1/omega/(rho_b*l*e+g2))))**(-1)
        return Q

    def ReynoldsNumber(self):
        # In 1962, Langlois derived the general form of Reynolds
        # equation based on Navier–Stokes equations and the general equations 
        # of viscous hydrodynamics. The Reynolds equation is obtained under 
        # the conditions that the modified Reynolds numbers Re are much smaller than unity.
        rho = self.fluid.density
        mu = self.fluid.dyna_viscosity
        width = self.beam.width
        omega = 2*np.pi*self.beam.freq # the radial frequency
        return rho*omega*width**2/4/mu


    def SaderCoefficient(self):
        # multiplicative correction coefficient to be applied to GammaFrontBack
        # in order to provide a more precise result in the case of infinitely 
        # thin rectangular beams [2] [Aoust's thesis page 87-88]
        # [2] it is accurate to within 0.1% over the entire range Re in [1E-6,1E4]
        Re = self.ReynoldsNumber()
        tau = np.log10(Re)
        SC1 = (0.91324-0.48274*tau+0.46842*tau**2-0.12886*tau**3+
            0.044055*tau**4-0.0035117*tau**5+0.00069085*tau**6)
        SC2 = (1-0.56964*tau+0.48690*tau**2-0.13444*tau**3+0.045155*tau**4-
            0.0035862*tau**5+0.00069085*tau**6)
        SC3 = (-0.024134-0.029256*tau+0.016294*tau**2-0.00010961*tau**3+
            0.000064577*tau**4-0.000044510*tau**5)
        SC4 = (1-0.59702*tau+0.55182*tau**2-0.18357*tau**3+0.079156*tau**4-
            0.014369*tau**5+0.0028361*tau**6)
        return SC1/SC2+1j*SC3/SC4


if __name__ == "__main__":
    Air = fm.Fluid(name='air')
    Air.info()
    Si = fm.Material(name='silicon_100')
    Si.info()

    Cant=be.Beam(material=Si,dico={'name':'cant',
                'length':150e-6,
                'width':22e-6,
                'thickness':4e-6,
                'gap':1.4e-6,
                'freq':160e3,
                'Qvac':10000})
    Cant.info()
    HyF = HydrodynamicFunction(Air,Cant)
    print('Re  : ',HyF.Re)
    print('Sad : ',HyF.Sader)
    print('Gfb : ',HyF.Gfb)
    print('Gup : ',HyF.Gup)
    print('Gsq : ',HyF.Gsqueeze)






