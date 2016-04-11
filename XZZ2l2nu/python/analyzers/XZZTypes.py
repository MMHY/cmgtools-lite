from CMGTools.XZZ2l2nu.analyzers.AutoFillTreeProducer  import * 



LLType = NTupleObjectType("LLType", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("mt",   lambda x : x.mt(), float),       
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
])



VVType = NTupleObjectType("VVType", baseObjectTypes=[], variables = [
  NTupleSubObject("LV",  lambda x : x['pair'],fourVectorType),
  NTupleVariable("deltaPhi",   lambda x : x['pair'].deltaPhi(), float),       
  NTupleVariable("deltaR",   lambda x : x['pair'].deltaR(), float),       
  NTupleVariable("mt",   lambda x : x['pair'].mt(), float),       
])


LLNuNuType = NTupleObjectType("LLNuNuType", baseObjectTypes=[VVType], variables = [
    NTupleSubObject("l1",  lambda x : x['pair'].leg1,LLType),
    NTupleSubObject("l1_l1",  lambda x : x['pair'].leg1.leg1,leptonTypeExtra),
    NTupleSubObject("l1_l2",  lambda x : x['pair'].leg1.leg2,leptonTypeExtra),
    NTupleSubObject("l2",  lambda x : x['pair'].leg2,metType),
])


JetType = NTupleObjectType("JetType", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("area",   lambda x : x.jetArea(), float),
    NTupleVariable("rawFactor",   lambda x : x.rawFactor(), float),
#    NTupleVariable("btag",   lambda x : x.bTag(), float),
#    NTupleVariable("nConstituents",   lambda x : len(x.constituents), int),
    # NTupleVariable("looseID",   lambda x : x.looseID, int),
    # NTupleVariable("tightID",   lambda x : x.tightID, int),
    # NTupleVariable("chargedHadronEnergyFraction",   lambda x : x.chargedHadronEnergyFraction(), float),
    # NTupleVariable("neutralHadronEnergyFraction",   lambda x : x.neutralHadronEnergyFraction(), float),
    # NTupleVariable("photonEnergyFraction",   lambda x : x.photonEnergyFraction(), float),
    # NTupleVariable("HFHadronEnergyFraction",   lambda x : x.HFHadronEnergyFraction(), float),
    # NTupleVariable("HFEMEnergyFraction",   lambda x : x.HFEMEnergyFraction(), float),
    # NTupleVariable("muonEnergyFraction",   lambda x : x.muonEnergyFraction(), float),
    # NTupleVariable("electronEnergyFraction",   lambda x : x.electronEnergyFraction(), float),
    # NTupleVariable("leptonEnergyFraction",   lambda x : x.leptonEnergyFraction(), float),

])


corrJetType = NTupleObjectType("corrJetType", baseObjectTypes=[JetType], variables = [
    NTupleVariable("jec_corr",   lambda x : x.corr, float), # JEC correction factor 
    NTupleVariable("jer_corr",   lambda x : x.corrJER if hasattr(x,"corrJER") else  1.0 ,float), # JER correction factor
])
