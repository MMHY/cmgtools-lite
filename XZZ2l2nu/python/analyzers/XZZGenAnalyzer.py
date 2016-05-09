from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsutils.genutils import *
from PhysicsTools.HeppyCore.statistics.counter import Counter, Counters
from CMGTools.XZZ2l2nu.tools.Pair import Pair
import math

class XZZGenAnalyzer( Analyzer ):
    """ Only select X->ZZ->2l2nu events

       """

    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(XZZGenAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

        self.filter = getattr(self.cfg_ana, 'filter', "None")
 
    def declareHandles(self):
        super(XZZGenAnalyzer, self).declareHandles()
        self.mchandles['genParticles'] = AutoHandle( 'prunedGenParticles', 'std::vector<reco::GenParticle>' )
                
    def beginLoop(self,setup):
        super(XZZGenAnalyzer,self).beginLoop(setup)

        self.counters.addCounter('XZZGenReport')
        self.count = self.counters.counter('XZZGenReport')
        self.count.register('XZZ2l2nu Events')
        self.count.register('XZZ2l2j Events')
        self.count.register('emu2nu Events')
        self.count.register('pass events')
        

    def makeMCInfo(self, event, bosonID = 23):
        verbose = getattr(self.cfg_ana, 'verbose', False)
        rawGenParticles = self.mchandles['genParticles'].product() 

        selectedGenParticles = []
        event.genZBosons = []
        event.genLeptons = []
        event.genElectrons = []
        event.genMuons = []
        event.genNeutrinos = []
        event.genJets = []
        event.genXZZ = []

        event.genZBosons = [ p for p in rawGenParticles if (abs(p.pdgId()) == bosonID) and p.numberOfDaughters() > 0 and abs(p.daughter(0).pdgId()) != 23 ]
        
        for zboson in event.genZBosons:
            selectedGenParticles.append(zboson)
            for i in xrange( zboson.numberOfDaughters() ):
                dau = zboson.daughter(i)
                dauid = dau.pdgId()
                if abs(dauid) in [11,13]:
                    event.genLeptons.append(dau)
                    selectedGenParticles.append(dau)
                    if abs(dauid)==11: event.genElectrons.append(dau)
                    if abs(dauid)==13: event.genMuons.append(dau)
                elif abs(dauid) in [12,14,16]:
                    event.genNeutrinos.append(dau)
                    selectedGenParticles.append(dau)
                elif abs(dauid) in range(7):
                    event.genJets.append(dau)
 
        event.genParticles = selectedGenParticles    
        
        if len(event.genZBosons)>=2 and len(event.genLeptons)>=2 and len(event.genNeutrinos)>=2:
            event.genIsXZZ2l2nu = True
        else: 
            event.genIsXZZ2l2nu = False
 
        if len(event.genZBosons)>=2 and len(event.genLeptons)>=2 and len(event.genJets)>=2:
            event.genIsXZZ2l2j = True
        else:
            event.genIsXZZ2l2j = False
            
        if len(event.genZBosons)>=2 and len(event.genElectrons)>=1 and len(event.genMuons)>=1 and len(event.genNeutrinos)>=2:
            event.genIsXZZemu2nu = True
        else:
            event.genIsXZZemu2nu = False

        if event.genIsXZZemu2nu:
            self.count.inc('emu2nu Events')

        if event.genIsXZZ2l2nu:
            self.count.inc('XZZ2l2nu Events') 

        if event.genIsXZZ2l2j:
            self.count.inc('XZZ2l2j Events')

        
        if len(event.genZBosons)>=2 :
            #for idx,genz in enumerate(event.genZBosons):
                #if genz.numberOfMothers()>1: 
                #    print 'GenZ ['+str(idx)+']: '
                #    print '   numberOfMothers(): '+str(genz.numberOfMothers())
            #    print 'GenZ ['+str(idx)+']: mother(0) '
            #    print genz.mother(0)

            #if event.genZBosons[0].mother(0) == event.genZBosons[1].mother(0):
            # if come from the same mother, take the mother
            #    event.genXZZ.append(event.genZBosons[0].mother(0))
            # not use the mother directly, because it doesn't seem to give the right mt... not yet understand why.     
            # create a new one sum the two Zs instead.
            #else:
            # if not, create a new one sum the two Zs
            genX = Pair(event.genZBosons[0],event.genZBosons[1],0)
            #print genX
            #print 'et='+str(genX.et())+'; Et='+str(genX.energy()*genX.pt()/genX.p())+'; mt='+str(genX.mt())+'; Mt='+str(math.sqrt((genX.energy()*genX.pt()/genX.p())**2-genX.pt()**2))
            #print 'pt='+str(genX.pt())+'; p='+str(genX.p())+'; E='+str(genX.energy())
            event.genXZZ.append({'pair':genX})
            

        if self.cfg_ana.verbose and event.genIsXZZ2l2nu:
            print "N gen Z bosons: "+str(len(event.genZBosons))
            print "N gen Leptons: "+str(len(event.genLeptons))
            print "N gen neutrinos: "+str(len(event.genNeutrinos))

    def process(self, event):
        self.readCollections( event.input )

        # if not MC, nothing to do
        if not self.cfg_comp.isMC: 
            return True
        # do MC level analysis
        self.makeMCInfo(event)
        
        # no gen filter
        if self.filter == "None":
            self.count.inc('pass events')
            return True
        # filter
        if self.filter=="ZZ2mu":
            if len(event.genMuons)>=2: 
                self.count.inc('pass events')
                return True
        elif self.filter=="ZZ2el":
            if len(event.genElectrons)>=2:
                self.count.inc('pass events')
                return True
        elif self.filter=="ZZ2l2nu":
            if event.genIsXZZ2l2nu:
                self.count.inc('pass events')
                return True

        # false if use filter but not pass any
        return False


import PhysicsTools.HeppyCore.framework.config as cfg
setattr(XZZGenAnalyzer,"defaultConfig",
    cfg.Analyzer(XZZGenAnalyzer,
        # Print out debug information
        verbose = False,
    )
)
