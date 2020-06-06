# Code to otimize a cantilever for photo-acoustic
# 2020-03-05 - First version
# 2020-03-06 - Import fluid and material from an other directory
# 2020-03-22 - Add Q support
# 2020-03-22 - Add a condition for mass_eff and hydrofunction caluculation
# 2020-03-23 - Add Q thermoelastic lcalculation


import numpy as np
from scipy.integrate import quad  # for integration
import sys
import os
# get current working directory
path = os.path.abspath(os.getcwd())
index = path[::-1].find('/')
path = path[:len(path)-index-1]
#print('current working directory : ',path)
# set the path of the module I want
sys.path.append(path) 

from hydrodynamic_function import hydro_function as hf
from fluid_material import fluid_material as fm
from geometry import beam as be

class Cantilever(be.Beam):
    def __init__(self,material,fluid,dico):
        super().__init__(material=material,dico=dico)
        self.material = material
        self.fluid = fluid
        if (self.length is not None) and (self.thickness is not None) and (self.width is not None):
            self.mass_eff = self.effective_mass()
        if dico['freq'] is None:
            self.freq = self.f0()
        # quality factor
        self.hydro = None
        self.Qsup = None # Q support
        self.Qvis = None # Q viscous
        self.Qted = None # Q thermoelastic
        self.Qaco = None # Q acoustic
        self.Qtot = None # Q total
        # suceptibility
        self.Xi = None # [m/N] or [s2/Kg]

    def effective_mass(self):
        return 0.248*self.mass

    def f0(self): 
        # resoance frequency in Hz for the first mode
        # for a clamped–free cantilever
        # red note-book page 2
        w = self.width
        t = self.thickness
        l = self.length
        m = self.mass
        I = w*t**3/12 # moment of inertia
        E = self.material.young # Material young's modulus
        mu = 1.875 #for the first harmonic
        freq = mu**2*np.sqrt(E*I/m/l**3)*1/(2*np.pi)
        return freq

    def length_unkonw(self):
        # deduce the length of the cnatilever from its frequency
        w = self.width
        t = self.thickness
        freq =self.freq # frequency in [Hz]
        I = (w*t**3)/12 # moment of inertia
        E = self.material.young # Material young's modulus
        rho = self.material.density # [kg/m3]
        mu = 1.875 #for the first harmonic
        l = (E*I/rho/w/t*(mu**2/2/np.pi/freq)**2)**(1./4)
        self.length = l
        self.mass = self.masse()
        self.mass_eff = self.effective_mass()   

    def phi(self,x,mode=1):
        # mode shape of a cantilever
        # mode=1 is for the first harmonic (ie the fundamental mode)
        alpha = [1.875,4.694,2.5*np.pi,3.5*np.pi,4.5*np.pi,5.5*np.pi,
                6.5*np.pi,7.5*np.pi,8.5*np.pi,9.5*np.pi][mode-1] 
        L = self.length
        X = alpha*x/L
        A = np.cosh(X)-np.cos(X)
        B = (np.sinh(alpha)-np.sin(alpha))/(np.cosh(alpha)+np.cos(alpha))
        C = np.sinh(X)-np.sin(X)
        return A-B*C 

    def Q_acoustic(self):
        # compute Q acoustic from Aoust thesis page 97
        # for code cf also red note book page 99
        e = self.thickness
        l = self.width
        L = self.length
        omega = 2*np.pi*self.freq
        c = self.fluid.speed_of_sound
        k = omega/c
        rho_p = self.material.density # [kg/m3]
        rho_f = self.fluid.density # [kg/m3]
        def Afunc(x):
            return self.phi(x,mode=1)**2
        A, _ = quad(Afunc, 0, L)
        A = e*A
        def Bfunc(phi):
            B1func = (np.sin(phi))**3 
            def B2func(phi):
                def B2Afunc(x):
                    return self.phi(x)*np.exp(-1j*k*x*np.cos(phi))
                temp, _ = np.abs(quad(B2Afunc, 0, L))**2
                return temp   
            return B1func*B2func(phi)
        B, _ = quad(Bfunc, 0, np.pi)
        Q = 256*rho_p/np.pi/rho_f/(k*l)**3*A/B
        self.Qaco = Q
        return Q

    def Q_info(self):
        print('QUALITY FACTOR :')
        print('Q support        = ',self.Qsup)
        print('Q viscous        = ',self.Qvis)
        print('Q vacuum         = ',self.Qvac)
        print('Q thermoelastic  = ',self.Qted)
        print('Q acoustic       = ',self.Qaco)
        print('Q TOTAL          = ',self.Qtot)

    def Q_support(self):
        # Formula from: 
        # Hao, Z., Erbil, A., & Ayazi, F. (2003). 
        # An analytical model for support loss in micromachined beam resonators 
        # with in-plane flexural vibrations. 
        # Sensors and Actuators, A: Physical, 109(1–2), 156–164. 
        # https://doi.org/10.1016/j.sna.2003.09.037
        beta = 0.597 
        # beta parameter for the first harmonic 
        #only of cantilever clamped free (cf red note book p94)
        nu = self.material.poisson
        L = self.length
        b = self.thickness
        Xi = (np.sin(np.pi*beta)-np.sinh(np.pi*beta))/(np.cos(np.pi*beta)+np.cosh(np.pi*beta))
        Q = 0.24*(1-nu)/(1+nu)/0.336/beta**2/Xi**2*(L/b)**3
        self.Qsup = Q 
        return Q  

    def Q_thermoelastic(self):
        # Q thermoelastic for a cantilever
        # Lifshitz, R., & Roukes, M. (2000). 
        # Thermoelastic damping in micro- and nanomechanical systems. 
        # Physical Review B - Condensed Matter and Materials Physics, 61(8), 5600–5609. 
        # https://doi.org/10.1103/PhysRevB.61.5600
        T = self.fluid.temperature
        rho = self.material.density # [kg/m3]
        E = self.material.young  # Young's modulus [N/m2]
        alpha = self.material.thermal_expansion # thermal expansion coefficient (K)
        Cp = self.material.Cp # heat capacity per unit volume at constant pressure [J/K/kg]
        K = self.material.thermal_conductivity # thermal conductivity [W/m/K]
        b = self.thickness 
        omega = 2*np.pi*self.freq
        ksi = b*np.sqrt(omega*rho*Cp/2/K)
        A = E*alpha**2*T/Cp/rho
        B = 6/ksi**2
        C = 6/ksi**3*(np.sinh(ksi)+np.sin(ksi))/(np.cosh(ksi)+np.cos(ksi))
        D = A*(B-C)
        Q = 1/D
        self.Qted = Q
        return Q

    def Q_viscous_sq(self):
        if self.hydro is None :         
            if (self.width is not None) and (self.thickness is not None):
               self.hydro = hf.HydrodynamicFunction(self.fluid,self)
            else:
                raise ValueError('Beam thickess and/or width should have a value.')
        # Q viscous with squeeze film
        self.hydro.Re = self.hydro.ReynoldsNumber()
        self.hydro.Sader = self.hydro.SaderCoefficient()
        self.hydro.Gfb = self.hydro.GammaFrontBack()
        self.hydro.Gup = self.hydro.GammaUpDown()
        self.hydro.Gsqueeze = self.hydro.GammaSqueeze()
        self.hydro.Ghydro = self.hydro.Sader*self.hydro.Gfb + self.hydro.Gup + self.hydro.Gsqueeze
        self.hydro.Q = self.hydro.QualityFactor()
        self.Qvis = self.hydro.Q
        return self.hydro.Q

    def Q_viscous_nosq(self):
        if self.hydro is None : 
            if (self.width is not None) and (self.thickness is not None):
                self.hydro = hf.HydrodynamicFunction(self.fluid,self)
            else:
                raise ValueError('Beam thickess and/or width should have a value.')
        # Q viscous without squeeze film
        self.hydro.Re = self.hydro.ReynoldsNumber()
        self.hydro.Sader = self.hydro.SaderCoefficient()
        self.hydro.Gfb = self.hydro.GammaFrontBack()
        self.hydro.Gup = self.hydro.GammaUpDown()
        self.hydro.Ghydro = self.hydro.Sader*self.hydro.Gfb + self.hydro.Gup 
        self.hydro.Q = self.hydro.QualityFactor()
        self.Qvis = self.hydro.Q
        return self.hydro.Q


    def surface(self):
        w = self.width
        l = self.length
        return l*w # surface in [m2]

    def Xi0(self):
        #suceptibility at resonance [m/N] or [s2/Kg]
        Q = self.Qtot
        f0 = self.freq
        m_eff = self.mass_eff
        Xi = Q/(2*np.pi*f0)**2/m_eff
        self.Xi = Xi
        return Xi


if __name__ == "__main__":
    Air = fm.Fluid(name='air')
    Air.info()
    print('speed of sound = %f m/s'%Air.speed_of_sound)
    Si = fm.Material(name='silicon_100')
    Si.info()

    Cant = Cantilever(material=Si,fluid=Air,dico={'name':'cant',
                'length':38.4e-3,
                'width':8e-3,
                'thickness':2e-3,
                'gap':None,
                'freq':4.76e3,
                'Qvac':None})
    Cant.Q_viscous_nosq()
    Cant.Q_support()
    Cant.Q_acoustic()
    Cant.info()
    Cant.Q_info()


