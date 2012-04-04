import FWCore.ParameterSet.Config as cms

electronHLTMatching = cms.EDProducer("ElectronHLTMatching",
                                     InputCollection = cms.InputTag("gsfElectrons"),
                                     TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                     HLTTriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                     TriggerPaths = cms.vstring("HLT_Ele8_v8",),
                                     RecoCuts = cms.string(""),
                                     HLTCuts = cms.string(""),
                                     DeltaR = cms.double(0.2),
                                     DoMatching = cms.bool(True)
                                     )

