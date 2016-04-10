import os
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.Heppy.analyzers.core.all import *
from PhysicsTools.Heppy.analyzers.objects.all import *
from PhysicsTools.Heppy.analyzers.gen.all import *
from PhysicsTools.HeppyCore.utils.deltar import *
from CMGTools.XZZ2l2nu.analyzers.Skimmer import *
from CMGTools.XZZ2l2nu.analyzers.XZZLeptonicVMaker import *
from CMGTools.XZZ2l2nu.analyzers.PackedCandidateLoader import *
from CMGTools.XZZ2l2nu.analyzers.XZZMultiFinalState  import *
from CMGTools.XZZ2l2nu.tools.leptonID  import *
from CMGTools.XZZ2l2nu.analyzers.XZZGenAnalyzer import *
from CMGTools.XZZ2l2nu.analyzers.XZZLeptonAnalyzer import *
from CMGTools.XZZ2l2nu.analyzers.XZZTriggerBitFilter import *
from CMGTools.XZZ2l2nu.analyzers.XZZVertexAnalyzer import *
from CMGTools.XZZ2l2nu.analyzers.XZZMETAnalyzer import *
from CMGTools.XZZ2l2nu.analyzers.XZZDumpEvtList import *
from CMGTools.XZZ2l2nu.analyzers.XZZJetAnalyzer import *

###########################
# define analyzers
###########################

skimAnalyzer = cfg.Analyzer(
    SkimAnalyzerCount, name='skimAnalyzerCount',
    useLumiBlocks = False,
    )

# Apply json file (if the dataset has one)
jsonAna = cfg.Analyzer(
    JSONAnalyzer, name="JSONAnalyzer",
    )

# Filter using the 'triggers' and 'vetoTriggers' specified in the dataset
triggerAna = cfg.Analyzer(
    XZZTriggerBitFilter, name="TriggerBitFilter",
    )

# This analyzer actually does the pile-up reweighting (generic)
pileUpAna = cfg.Analyzer(
    PileUpAnalyzer, name="PileUpAnalyzer",
    true = True,  # use number of true interactions for reweighting
    makeHists=False
    )


## Gen Info Analyzer 
genAna = cfg.Analyzer(
    XZZGenAnalyzer, name="XZZGenAnalyzer",
    # Print out debug information
    verbose = False,
    filter = "None",
    )

# Select a list of good primary vertices (generic)
vertexAna = cfg.Analyzer(
    XZZVertexAnalyzer, name="VertexAnalyzer",
    allVertices = "offlineSlimmedPrimaryVertices",
    vertexWeight = None,
    fixedWeight = 1,
    verbose = False
    )

lepAna = cfg.Analyzer(
    XZZLeptonAnalyzer, name="leptonAnalyzer",
    muons='slimmedMuons',
    electrons='slimmedElectrons',
    packedCandidates = 'packedPFCandidates',
    rhoMuon= 'fixedGridRhoFastjetCentralNeutral',
    rhoElectronMiniIso = 'fixedGridRhoFastjetCentralNeutral',
    rhoElectronPfIso = 'fixedGridRhoFastjetAll',
    applyIso = True,
    electronIDVersion = 'looseID', # can be looseID or HEEPv6
    electronIsoVersion = 'pfISO', # can be pfISO or miniISO
    mu_isoCorr = "rhoArea" ,
    ele_isoCorr = "rhoArea" ,
    mu_effectiveAreas = "Spring15_25ns_v1",
    ele_effectiveAreas = "Spring15_25ns_v1",
    miniIsolationPUCorr = None, # Allowed options: 'rhoArea' (EAs for 03 cone scaled by R^2), 'deltaBeta', 
                                     # 'raw' (uncorrected), 'weights' (delta beta weights; not validated)
                                     # Choose None to just use the individual object's PU correction
    )

## Jets Analyzer (generic)
jetAna = cfg.Analyzer(
    JetAnalyzer, name='jetAnalyzer',
    debug=False,
    jetCol = 'slimmedJets',
    copyJetsByValue = False,      #Whether or not to copy the input jets or to work with references (should be 'True' if JetAnalyzer is run more than once)
    genJetCol = 'slimmedGenJets',
    rho = ('fixedGridRhoFastjetAll','fixedGridRhoAll',''), # it was ('fixedGridRhoFastjetAll','','') 
    jetPt = 5., # default used to be 25.
    jetEta = 4.7,
    jetEtaCentral = 2.4,
    jetLepDR = 0.4,
    jetLepArbitration = (lambda jet,lepton : lepton), # you can decide which to keep in case of overlaps; e.g. if the jet is b-tagged you might want to keep the jet
    cleanSelectedLeptons = True, #Whether to clean 'selectedLeptons' after disambiguation. Treat with care (= 'False') if running Jetanalyzer more than once
    minLepPt = 10,
    relaxJetId = False,
    doPuId = False, # Not commissioned in 7.0.X
    recalibrateJets = True, #'MC', # True, False, 'MC', 'Data'
    applyL2L3Residual = True, # Switch to 'Data' when they will become available for Data
    recalibrationType = "AK4PFchs",
    mcGT     = "Fall15_25nsV2_MC", # txt file pattern used in 74X: 'Summer15_25nsV6_MC'
    dataGT   = "Fall15_25nsV2_DATA", #Summer15_25nsV6_DATA
    mcGT_jer    = "Summer15_25nsV6_MC", 
    dataGT_jer  = "Summer15_25nsV6_DATA", 
    jecPath = "${CMSSW_BASE}/src/CMGTools/XZZ2l2nu/data/jec/", #${CMSSW_BASE}/src/CMGTools/RootTools/data/jec/
    jerPath = "${CMSSW_BASE}/src/CMGTools/XZZ2l2nu/data/jer/", #${CMSSW_BASE}/src/CMGTools/RootTools/data/jec/
    shiftJEC = 0, # set to +1 or -1 to apply +/-1 sigma shift to the nominal jet energies
    addJECShifts = False, # if true, add  "corr", "corrJECUp", and "corrJECDown" for each jet (requires uncertainties to be available!)
    smearJets = False,
    shiftJER = 0, # set to +1 or -1 to get +/-1 sigma shifts  
    alwaysCleanPhotons = False,
    cleanGenJetsFromPhoton = False,
    cleanJetsFromFirstPhoton = False,
    cleanJetsFromTaus = False,
    cleanJetsFromIsoTracks = False,
    doQG = False,
    do_mc_match = False,
    collectionPostFix = "",
    calculateSeparateCorrections = True, # should be True if recalibrateJets is True, otherwise L1s will be inconsistent
    calculateType1METCorrection  = True,
    type1METParams = { 'jetPtThreshold':15., 'skipEMfractionThreshold':0.9, 'skipMuons':True }, # numbers for AK4CHS jets
    )

metAna = cfg.Analyzer(
    XZZMETAnalyzer, name="metAnalyzer",
    metCollection     = "slimmedMETs",
    noPUMetCollection = "slimmedMETs",
    copyMETsByValue = False,
    doTkMet = False,
    doMetNoPU = True,
    doMetNoMu = False,
    doMetNoEle = False,
    doMetNoPhoton = False,
    recalibrate = "type1", # or "type1", or True, or False
    applyJetSmearing = True, # does nothing unless the jet smearing is turned on in the jet analyzer
    old74XMiniAODs = False, # set to True to get the correct Raw MET when running on old 74X MiniAODs
    jetAnalyzerPostFix = "",
    candidates='packedPFCandidates',
    candidatesTypes='std::vector<pat::PackedCandidate>',
    dzMax = 0.1,
    collectionPostFix = "",
    )


leptonicVAna = cfg.Analyzer(
    XZZLeptonicVMaker,
    name='leptonicVMaker',
    selectMuMuPair = (lambda x: (x.leg1.highPtID or x.leg2.highPtID) and ((x.leg1.pt()>50.0 and abs(x.leg1.eta())<2.1) or (x.leg2.pt()>50.0 and abs(x.leg2.eta())<2.1))),
    selectElElPair = (lambda x: x.leg1.pt()>115.0 or x.leg2.pt()>115.0 ),
    selectVBoson = (lambda x: x.pt()>50.0 and x.mass()>60.0 and x.mass()<120.0) 
    )

packedAna = cfg.Analyzer(
    PackedCandidateLoader,
    name = 'PackedCandidateLoader',
    select=lambda x: x.pt()<13000.0
    )

multiStateAna = cfg.Analyzer(
    XZZMultiFinalState,
    name='MultiFinalStateMaker',
    selectPairLLNuNu = (lambda x: x.leg1.pt()>50.0 and x.leg1.mass()>60.0 and x.leg1.mass()<120.0 and x.leg2.pt()>50.0),
    suffix = '',
    )

# Create flags for MET filter bits

eventFlagsAna = cfg.Analyzer(
    TriggerBitAnalyzer, name="EventFlags",
    processName = 'PAT',
    fallbackProcessName = 'RECO',
    outprefix   = 'Flag',
    triggerBits = {
        "HBHENoiseFilter" : [ "Flag_HBHENoiseFilter" ],
        "HBHENoiseIsoFilter" : [ "Flag_HBHENoiseIsoFilter" ],
        "CSCTightHaloFilter" : [ "Flag_CSCTightHaloFilter" ],
        "hcalLaserEventFilter" : [ "Flag_hcalLaserEventFilter" ],
        "EcalDeadCellTriggerPrimitiveFilter" : [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ],
        "goodVertices" : [ "Flag_goodVertices" ],
        "trackingFailureFilter" : [ "Flag_trackingFailureFilter" ],
        "eeBadScFilter" : [ "Flag_eeBadScFilter" ],
        "ecalLaserCorrFilter" : [ "Flag_ecalLaserCorrFilter" ],
        "trkPOGFilters" : [ "Flag_trkPOGFilters" ],
        "trkPOG_manystripclus53X" : [ "Flag_trkPOG_manystripclus53X" ],
        "trkPOG_toomanystripclus53X" : [ "Flag_trkPOG_toomanystripclus53X" ],
        "trkPOG_logErrorTooManyClusters" : [ "Flag_trkPOG_logErrorTooManyClusters" ],
        "METFilters" : [ "Flag_METFilters" ],
    }
    )

# Create flags for trigger bits
triggerFlagsAna = cfg.Analyzer(
    TriggerBitAnalyzer, name="TriggerFlags",
    processName = 'HLT',
    triggerBits = {
    }
    )

dumpEvents = cfg.Analyzer(
    XZZDumpEvtList, name="XZZDumpEvtList",
    )


################
coreSequence = [
    skimAnalyzer,
    genAna,
    jsonAna,
    triggerAna,
    pileUpAna,
    vertexAna,
    lepAna,
    jetAna,
    metAna,
    leptonicVAna,
#    packedAna,
    multiStateAna,
    eventFlagsAna,
    triggerFlagsAna
]

###################




