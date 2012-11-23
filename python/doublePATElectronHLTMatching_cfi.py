import FWCore.ParameterSet.Config as cms

doublePATElectronHLTMatching = cms.EDProducer("DoublePATElectronHLTMatching",
                                           InputCollection = cms.InputTag("patElectronsWithTrigger"),
                                           TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                           HLTTriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                           TriggerPaths = cms.vstring(
                                                                      #"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
                                                                      "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*",
                                                                      #"HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*",
                                                                      ),
                                           RecoCuts   = cms.string(""),
                                           HLTCuts    = cms.string(""),
                                           TagLeg     = cms.bool(True),
                                           DeltaR     = cms.double(0.2),
                                           DoMatching = cms.bool(True)
                                           )

