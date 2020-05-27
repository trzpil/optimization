# Give physical parameter for different fluid and material
# 2020-03-06 - First version
# 2020-03-20 - Add Water
# 2020-04-09 - Add for Air speed of sound and gamma
# 2020-05-27 - Add gas species

import numpy as np

class GasSpecie():
    # contain the physical parameters of the target gas
    def __init__(self,dico):
        self.name = dico['name']
        self.wavelength = dico['abs wavelength [m]']
        self.pressure = dico['pressure [Pa]']
        self.temperature = dico['temperature [K]']
        self.cross_section = dico['cross section [cm2/molecule]']
        self.relaxation_time = dico['relaxation time [s]']
        self.concentration = dico['concentration [Nbr_mlc/Nbr_tot]']
        self.comp_gas_density()
        self.comp_absorption_coeff()
        self.comp_absorption()

    def comp_absorption(self):
        # compute absorption in % 
        # for a optical path of 1m
        length = 100 # [cm]
        self.absorption = 1-np.exp(-self.absorption_coeff*length)
        self.absorbance = self.absorption_coeff*length

    def comp_absorption_coeff(self):
        # compute the absorption coefficient [2]
        self.absorption_coeff =  self.cross_section * self.gas_densty 

    def comp_gas_density(self):
        # see ref [2] Loschmidt number
        N_L = 2.479e19 # [mol/cm3/atm]
        # compute N_tot the total density of molecules 
        # at temperature T and pressure P
        P = self.pressure/101325 # in atmosphere
        T = self.temperature
        N_tot = N_L*296/T*P
        # compute N the gas density in [mol/cm3]
        C = self.concentration
        N = C*N_tot
        self.gas_densty = N

    def info(self):
        print('GAS SPECIE '+ self.name+' : ')
        print('absorption wavelength = %f um'
            %(self.wavelength*1e6))
        print('pressure = %g Pa'
            %(self.pressure))
        print('temperature = %.1f K'
            %(self.temperature))
        print('cross section = %g cm2/molecule'
            %(self.cross_section))
        print('relaxation time = %f us'
            %(self.relaxation_time*1e6))
        print('concentration = %g ppmv'
            %(self.concentration*1e6))
        print('gas density = %g mol/cm3'
            %(self.gas_densty))
        print('absorption coefficient = %g cm-1'
            %(self.absorption_coeff))
        print('absorption (on 1m) = {:.0%} '
            .format(self.absorption))
        print('absorbance (on 1m) = {:.0%} '
            .format(self.absorbance))

class Fluid():
    def __init__(self,name=None):
        ## INTRINSIC
        self.density = None # fluid density [kg/m3]
        self.dyna_viscosity = None # mu dynamic viscosity [kg/(m.s)]
        ## EXTRINSIC
        self.pressure = 101325 # Sea level standard atmospheric pressure [Pa]
        self.temperature = 300 # in [K]
        ## PREDEFINED FLUID
        self.name = name
        if self.name in 'air':
            self.fluid_is_air()
        if self.name in 'water':
            self.fluid_is_water()

    def fluid_is_air(self):
        def compute_density(Rspe = 287.058):
            # Rpse The specific gas constant for dry air is 287.058 [J/(kg·K)] SI units
            # P the pressure in [Pa]
            P = self.pressure
            # T temperature in [K]
            T = self.temperature 
            return P/T/Rspe #[kg/m3]
        def compute_dyna_visco():
            # T temperature in [K]
            T = self.temperature
            # mu dynamic viscosity [kg/(m.s)] https://fr.wikipedia.org/wiki/Air
            mu = 8.8848e-15*T**3-3.2398e-11*T**2+6.2657e-8*T+2.3543e-6
            return mu # [kg/(m.s)]
        self.density = compute_density()
        self.dyna_viscosity = compute_dyna_visco()
        self.speed_of_sound = 20.05*np.sqrt(self.temperature) # [m/s]
        # /!\ ONLY AT 300K 
        self.relative_permittivity = 1.000536 # relative permittivity [F/m] = [C2/(N m2)]
        self.Cp = 1005 #Isobaric specific heat  is used for air in a constant pressure [J/Kg/K]
        self.Cv = 718 #Isochoric specific heat is used for air in a constant-volume [J/Kg/K]
        self.Gamma = self.Cp/self.Cv  #Heat capacity ratio, also known as the adiabatic index, the ratio of specific heats,

    def fluid_is_water(self):
        self.density = 1000 #[kg/m3]
        self.dyna_viscosity = 0.0010#518 # @18°C mu dynamic viscosity [kg/(m.s)]

    def info(self):
        print('FLUID '+ self.name+' : ')
        print('temperature = %g K'%(self.temperature))
        print('pressure = %g Pa'%(self.pressure))
        print('density = %g kg/m3 at %g K and %g Pa '
            %(self.density,self.temperature,self.pressure))
        print('dynamic viscosity = %g kg/m/s at 300 K '
            %(self.dyna_viscosity))
        if hasattr(self, 'Gamma'): 
            print('air heat capacity ratio : ',self.Gamma)
        if hasattr(self, 'speed_of_sound'): 
            print('speed of sound : %f m/s'%self.speed_of_sound)

class Material():
    def __init__(self,name=None):
        # PREDEFINED MATERIAL
        self.name = name
        if self.name in 'silicon':
            self.material_is_silicon()
        if self.name in 'silicon_100':
            self.material_is_silicon_100()
        if self.name in 'silicon_110':
            self.material_is_silicon_110()        
        if self.name in 'quartz':
            self.material_is_quartz()
        if self.name in 'diamond':
            self.material_is_diamond()
        

    def material_is_diamond(self):
        self.density = 3530 # [kg/m3]

    def material_is_quartz(self):
        # density for alpha-quartz
        self.density = 2648 # [kg/m3]

    def material_is_silicon(self):
        self.density = 2330 # [kg/m3]

    def material_is_silicon_100(self):
        #<100> directions ("45° off-axis")
        #For the “off-axis” direction (45◦ diagonal to flat), use E100 = 130 GPa [1]
        self.density = 2330 # [kg/m3]
        self.young = 130e9 # Young's modulus [N/m2]
        self.poisson = 0.28 # Poisson's ratio
        self.thermal_expansion = 2.6e-6 # thermal expansion coefficient (1/K)
        self.Cp = 700 # heat capacity per unit volume at constant pressure [J/K/kg]
        self.thermal_conductivity = 90  # thermal conductivity [W/m/K)]

    def material_is_silicon_110(self):
        #<110> directions ("X or Y axis")
        #For the “X-or Y-axis” direction (parallel/perpendicular to flat), use E110 = 169 GPa [1]      
        self.density = 2330 # [kg/m3]
        self.young = 169e9 # Young's modulus [N/m2]           

    def info(self):
        print('MATERIAL '+ self.name+' : ')
        print('density = %g kg/m3'
            %(self.density))
        print('Young’s modulus = %g GPa'
            %(self.young*1e-9))

if __name__ == "__main__":
    Air = Fluid(name='air')
    Air.info()

    Ch4 = GasSpecie(dico={
        'name':'CH4',
        'abs wavelength [m]':1.65e-6,
        'pressure [Pa]':Air.pressure,
        'temperature [K]':Air.temperature,
        'cross section [cm2/molecule]':1.55e-20,
        'relaxation time [s]':15e-6,
        'concentration [Nbr_mlc/Nbr_tot]':0.01
        })
    Ch4.info()

    Si = Material(name='silicon_100')
    Si.info()
    # print(Si.poisson)


#################################################
###   REF
#################################################

# [1]
# Hopcroft, M. A., Nix, W. D., & Kenny, T. W. (2010). 
# "What is the Young’s modulus of silicon? "
# Journal of Microelectromechanical Systems, 19(2), 229–238. 
# https://doi.org/10.1109/JMEMS.2009.2039697  
# [2]
# Jean-Philippe BESSON
# "photoacoustic spectroscopy for multi-gas sensing using near infrared lasers"
# Thesis
