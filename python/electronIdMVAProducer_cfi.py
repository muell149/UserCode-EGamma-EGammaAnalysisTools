import FWCore.ParameterSet.Config as cms

basePath = '/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/'

mvaTrigV0 = cms.EDFilter("ElectronIdMVAProducer",
                            verbose = cms.untracked.bool(False),
                            vertexTag = cms.InputTag('offlinePrimaryVertices'),
                            electronTag = cms.InputTag('gsfElectrons'),
                            reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB'),
                            reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE'),
                            method = cms.string("BDTCat_BDTG_TrigV0"),
                            #mvaWeightFile = cms.string("EGamma/EGammaAnalysisTools/data/Electrons_BDTGCat_TrigV0.weights.xml"),
                            mvaWeightFile = cms.string(basePath + "Electrons_BDTGCat_TrigV0.weights.xml"),
                            Trig = cms.bool(True),
)

mvaNonTrigV0 = cms.EDFilter("ElectronIdMVAProducer",
                            verbose = cms.untracked.bool(False),
                            vertexTag = cms.InputTag('offlinePrimaryVertices'),
                            electronTag = cms.InputTag('gsfElectrons'),
                            reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB'),
                            reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE'),
                            method = cms.string("BDTCat_BDTG_NonTrigV0"),
                            #mvaWeightFile = cms.string("EGamma/EGammaAnalysisTools/data/Electrons_BDTGCat_NonTrigV0.weights.xml"),
                            mvaWeightFile = cms.string(basePath + "Electrons_BDTGCat_NonTrigV0.weights.xml"),
                            Trig = cms.bool(False),
)
