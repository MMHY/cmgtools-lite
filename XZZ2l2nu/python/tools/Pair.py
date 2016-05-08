import math

from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi


class Pair(object):
    def __init__(self,leg1,leg2,pdg = 0):
        self.leg1 = leg1
        self.leg2 = leg2
        self.pdg = pdg
        self.LV = leg1.p4()+leg2.p4()

        et1 = math.sqrt(leg1.mass()*leg1.mass()+leg1.pt()*leg1.pt())
        et2 = math.sqrt(leg2.mass()*leg2.mass()+leg2.pt()*leg2.pt())
        self.MT  =math.sqrt(self.leg1.p4().mass()*self.leg1.p4().mass()+\
            self.leg2.p4().mass()*self.leg2.p4().mass()+2*(et1*et2-self.leg1.p4().px()*self.leg2.p4().px()-self.leg1.p4().py()*self.leg2.p4().py()))

        et2a = math.sqrt(91.188**2+leg2.pt()*leg2.pt())  
        self.MTa = math.sqrt(self.leg1.p4().mass()*self.leg1.p4().mass()+\
            91.188**2+2*(et1*et2a-self.leg1.p4().px()*self.leg2.p4().px()-self.leg1.p4().py()*self.leg2.p4().py()))

        et1b = leg1.p4().Et()
        et2b = math.sqrt(91.188**2+leg2.pt()**2)
        self.MTb = math.sqrt( (et1b+et2b)**2 - (leg1.pt()**2+leg2.pt()**2+2.0*(self.leg1.p4().px()*self.leg2.p4().px()+self.leg1.p4().py()*self.leg2.p4().py())) ) 

        et1c = math.sqrt((leg1.mass())**2+leg1.pt()**2)
        et2c = math.sqrt((leg1.mass())**2+leg2.pt()**2)
        self.MTc = math.sqrt( (et1c+et2c)**2 - (leg1.pt()**2+leg2.pt()**2+2.0*(self.leg1.p4().px()*self.leg2.p4().px()+self.leg1.p4().py()*self.leg2.p4().py())) )

        et1d = leg1.energy()*leg1.pt()/leg1.p()
        et2d = leg2.energy()*leg2.pt()/leg2.p()
        self.MTd = math.sqrt( (et1d+et2d)**2 - self.LV.pt()**2 )

        #et1e = math.sqrt(91.188**2+leg1.pt()*leg1.pt())
        #et2e = math.sqrt(91.188**2+leg2.pt()*leg2.pt())
        #self.MTe = math.sqrt( (et1e+et2e)**2 - self.LV.pt()**2 )

    def rawP4(self):
        return self.leg1.p4()+self.leg2.p4()

    def p4(self):
        return self.LV
    
    def m(self):
        return self.LV.mass()

    def e(self):
        return self.LV.energy()
    
    def pdgId(self):
        return self.pdg
    
    def mt2(self):
        return self.MT*self.MT

    def mt(self):
        return self.MT
  
    def mta(self):
        return self.MTa

    def mtb(self):
        return self.MTb

    def mtc(self):
        return self.MTc
 
    def mtd(self):
        return self.MTd

    #def mte(self):
    #    return self.MTe

    def pt(self):
        return self.LV.pt()

    def p(self):
        return self.LV.p()

    def px(self):
        return self.LV.px()

    def py(self):
        return self.LV.py()

    def pz(self):
        return self.LV.pz()

    def et(self):
        return self.LV.Et()

    def eta(self):
        return self.LV.eta()

    def phi(self):
        return self.LV.phi()

    def deltaPhi(self):
        return abs(deltaPhi(self.leg1.phi(),self.leg2.phi()))

    def deltaR(self):
        return abs(deltaR(self.leg1.eta(),self.leg1.phi(),self.leg2.eta(),self.leg2.phi()))

    def __getattr__(self, name):
        return getattr(self.LV,name)

