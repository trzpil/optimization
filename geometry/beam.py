# Code to otimize a cantilever for photo-acoustic
# 2020-03-05 - First version
# 2020-03-06 - Import fluid and material from an other directory
# 2020-03-22 - Add a condition for mass caluculation


class Beam():
    def __init__(self,material,dico):
        self.name = dico['name']
        ## MATERAIL
        self.density = material.density # [kg/m3]
        ## GEOMETRY
        self.length = dico['length'] # length [m]
        # Q formula only implicitly depends on the length of the beam 
        # through the resonant frequency
        self.width = dico['width'] # width [m]
        self.thickness = dico['thickness'] # thickness in [m]
        self.gap = dico['gap'] # gap of air for the squeeze film
        if (self.length is not None) and (self.thickness is not None)  and (self.width is not None):
            self.mass = self.masse()
        else:
            self.mass = None
        ## RESONANCE
        self.freq = dico['freq'] # [Hz]
        self.Qvac = dico['Qvac'] # quality factor under vacuum

    def info(self):
        print('BEAM '+ self.name+' : ')
        if self.length is not None:
            print('length = %g mm'%(self.length*1e3))
        else:
            print('length = None')
        if self.width is not None:
            print('width = %g um'%(self.width*1e6))
        else:
            print('width = None')
        if self.thickness is not None:
            print('thickness = %g um'%(self.thickness*1e6))
        else:
            print('thickness = None')
        if self.gap is not None:
            print('gap = %g um'%(self.gap*1e6))
        else:
            print('gap = None')
        if self.freq is not None:
            print('resonant frequency = %g kHz'%(self.freq/1e3))
        else:
            print('resonant frequency = None')
        if self.Qvac is not None:
            print('quality factor under vacum = %g'%(self.Qvac))
        else:
            print('quality factor under vacum = None')
        if self.mass is not None:
            print('mass = %g ug'%(self.mass*1e9))
        else:
            print('mass = None')
        print('\n')

    def masse(self):
        m = self.density*self.length*self.width*self.thickness #[kg]
        return m 