from CMGTools.XZZ2l2nu.analyzers.AutoFillTreeProducer  import * 
from CMGTools.XZZ2l2nu.analyzers.XZZTypes  import * 
from CMGTools.XZZ2l2nu.analyzers.Skimmer  import * 
import PhysicsTools.HeppyCore.framework.config as cfg

vvSkimmer = cfg.Analyzer(
    Skimmer,
    name='vvSkimmer',
    required = ['LLNuNu', 'ElMu']
)

leptonSkimmer = cfg.Analyzer(
    Skimmer,
    name='leptonSkimmer',
    required = ['inclusiveLeptons']
)

jetTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='jetTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"), 
     ],
     globalObjects =  {
         "met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, after full type 1 corrections"),
         "met_miniAod" : NTupleObject("met_raw", metType, help="PF E_{T}^{miss}, after type 1 corrections in miniAOD"),
         "met_JEC" : NTupleObject("met_JEC", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC added"),
         "met_JECUp" : NTupleObject("met_JECUp", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC up"),
         "met_JECDown" : NTupleObject("met_JECDown", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC down"),
         "met_JECJER" : NTupleObject("met_JECJER", metType, help="PF E_{T}^{miss}, after type 1 corrections with JEC+JER jets"),
     },

     collections = {
         "selectedLeptons" : NTupleCollection("lep",leptonType,10, help="selected leptons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         "genZBosons" : NTupleCollection("genZ", genParticleType, 10, help="Generated V bosons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         "genZBosons" : NTupleCollection("genZ", genParticleType, 10, help="Generated V bosons"),
         "jets"  : NTupleCollection("jet_corr",corrJetType,15, help="all jets with new JEC for 76X applied"),
         "jets_raw"  : NTupleCollection("jet",JetType,15, help="all jets from miniAOD"),
     }
)

vvTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='vvTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
         NTupleVariable("nLL",lambda ev: len(ev.LL) , int),       
         NTupleVariable("nElMu",lambda ev: len(ev.ElMu) , int),       
         NTupleVariable("nLLNuNu",lambda ev: len(ev.LLNuNu) , int),       
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"), 
         NTupleVariable("triggersf",  lambda x : getattr(x,'trgsf',1), help="singleelectron/muon trigger sf"),
         NTupleVariable("triggersfUp",  lambda x : getattr(x,'trgsfUp',1), help="singleelectron/muon trigger sf upper"),
         NTupleVariable("triggersfLo",  lambda x : getattr(x,'trgsfLo',1), help="singleelectron/muon trigger sf lower"),
     ],
     globalObjects =  {
         "met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, after type 1 corrections"),
         #"met_miniAod" : NTupleObject("met_raw", metType, help="PF E_{T}^{miss}, after type 1 corrections in miniAOD"),
         #"met_JEC" : NTupleObject("met_JEC", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC added"),
         "met_JECUp" : NTupleObject("met_JECUp", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC up"),
         "met_JECDown" : NTupleObject("met_JECDown", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC down"),
         #"met_JECJER" : NTupleObject("met_JECJER", metType, mcOnly=True, help="PF E_{T}^{miss}, after type 1 corrections with JEC+JER jets"),
     },

     collections = {
         #"LL"  : NTupleCollection("Zll",LLType,5, help="Z to ll"),
         "ElMu"  : NTupleCollection("elmu",LLType,5, help="electron - muon pair for non-resonant bkg"),
         "selectedLeptons" : NTupleCollection("lep",leptonType,10, help="selected leptons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         #"genZBosons" : NTupleCollection("genZ", genParticleType, 10, help="Generated V bosons"),
         "LLNuNu"     : NTupleCollection("llnunu",LLNuNuType ,5, help="VV candidate with di-lepton and MET"),
         #"genXZZ" : NTupleCollection("genX", VVType, 10, mcOnly=True, help="Generated X->ZZ"),
         "LHE_weights"    : NTupleCollection("LHEweight",  weightsInfoType, 1000, help="LHE weight info"),
         #"jets"       : NTupleCollection("jet_corr",corrJetType,15, help="all jets with new JEC for 76X applied (JER corrected if isMC)"),
         #"jets_raw"   : NTupleCollection("jet",JetType,15, help="all jets from miniAOD"),
     }
)

lepeffTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='lepeffTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),
     ],
     globalObjects =  {
         "met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, after type 1 corrections"),
     },
     collections = {
         "llpair"  : NTupleCollection("llpair",llpairType ,5, help="lepton eff study"),
     }
)

leptonTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='leptonTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
#     PDFWeights = PDFWeights,
#     globalVariables = susyMultilepton_globalVariables,
#     globalObjects = susyMultilepton_globalObjects,
     collections = {
            "genleps"          : NTupleCollection("gen",     genParticleWithLinksType, 10, help="Generated leptons (e/mu) from W/Z decays"),                                                                                                
            "inclusiveLeptons" : NTupleCollection("l",    leptonTypeExtra, 10, help="Inclusive Leptons"),                                                                                                
     }
)
