import ROOT
import random
import math
from  itertools import combinations
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import *
import PhysicsTools.HeppyCore.framework.config as cfg
from CMGTools.XZZ2l2nu.tools.Pair import *

class XZZMuonEffTree( Analyzer ):

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZMuonEffTree,self).__init__(cfg_ana,cfg_comp,looperName)
        self.genfilter=getattr(self.cfg_ana,"genfilter",False)
        self.eithercharge=getattr(self.cfg_ana,"eithercharge",False)
        self.checktag=getattr(cfg_ana,'checktag',False)
        self.pfbkg=getattr(cfg_ana,'pfbkg',False)
        self.pctpf=getattr(cfg_ana,'pctpf',{3:0,4:0})
        self.muHLT=getattr(cfg_ana,'muHLT','')
        self.npf={3:0,4:0,"all":.01}



    def declareHandles(self): 
        super(XZZMuonEffTree, self).declareHandles()
        if self.checktag:
            self.handles['trgresults']=AutoHandle(("TriggerResults","","HLT"),"edm::TriggerResults")
            self.handles['selectedtrg']=AutoHandle('selectedPatTrigger','std::vector<pat::TriggerObjectStandAlone>')
        if self.pfbkg:
            self.handles['pfcandidate']=AutoHandle("packedPFCandidates","std::vector<pat::PackedCandidate>")

    def checkgen(self,lep):
        try:
            mom=lep.physObj.genParticle()
            while mom.pdgId()!=23:
                mom=mom.mother()
            return True
        except:
            return False

    def istag(self,mu,tobs):
        itg=1
        if not mu.muonID("POG_ID_Tight"): itg-=1
#        if (mu.physObj.pfIsolationR04().sumChargedHadronPt + max(0., mu.physObj.pfIsolationR04().sumNeutralHadronEt + mu.physObj.pfIsolationR04().sumPhotonEt - 0.5*mu.physObj.pfIsolationR04().sumPUPt))/mu.pt()>.2: itg-=2
        if mu.absIsoWithFSR()/mu.pt()>.2: itg-=2
        if not [i for i in tobs if deltaR(mu.eta(),mu.phi(),i.eta(),i.phi())<.3]: itg-=4
        if mu.pt()<21: itg-=8
        return itg

    def process(self, event):
        self.readCollections( event.input )
        event.llpair=[]
        if self.pfbkg:
            event.selectedhadrons=[i for i in self.handles['pfcandidate'].product() if i.pt()>20 and abs(i.pdgId())==211 and abs(i.eta())<2.4]
            event.selectedhadrons.sort(key = lambda l : l.pt(), reverse = True)
            if len(event.selectedhadrons)>4 and self.npf[4]/self.npf['all']<self.pctpf[4]:
                event.zllcandidates=event.selectedhadrons[::len(event.selectedhadrons)/4][:4]
                self.npf[4]+=1.
            elif len(event.selectedhadrons)>3 and self.npf[3]/self.npf['all']<self.pctpf[3]:
                event.zllcandidates=event.selectedhadrons[::len(event.selectedhadrons)/3][:3]
                self.npf[3]+=1.
            elif len(event.selectedhadrons)>2:event.zllcandidates=event.selectedhadrons[::len(event.selectedhadrons)/2][:2]
            else: return False
            self.npf['all']+=1.
        else:
            event.zllcandidates=[i for i in event.selectedMuons if i.track().isNonnull()]
            if len(event.zllcandidates)<2: return False
        for l1,l2 in combinations(event.zllcandidates,2):
            if l1.pdgId() == -l2.pdgId() or self.eithercharge:
                pair = Pair(l1,l2,23)
                if abs(l1.dz()-l2.dz())<.4 or self.pfbkg: event.llpair.append(pair)
        if not event.llpair: return False
        event.llpair=[min(event.llpair,key = lambda x: abs(x.M()-91.118))]
        if self.checktag:
            names = event.input.object().triggerNames(self.handles['trgresults'].product())
            tobs=[]
            for i in self.handles['selectedtrg'].product():
                i.unpackPathNames(names)
                pNames=list(i.pathNames())
                if self.muHLT and [pN for pN in pNames if self.muHLT in pN]:tobs.append(i)
            for i in event.llpair:
                i.leg1.istag=self.istag(i.leg1,tobs)
                i.leg2.istag=self.istag(i.leg2,tobs)
        if self.genfilter and self.cfg_comp.isMC:
            for i in event.llpair:
                i.leg1.xdaughter=self.checkgen(i.leg1)
                i.leg2.xdaughter=self.checkgen(i.leg2)
        return True




        
            

        


                
                
