import FWCore.ParameterSet.Config as cms

process = cms.Process("ExISO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C61C71DC-8061-E111-AAEB-0025B32036D2.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C49AC4E5-7F61-E111-8FB2-003048F1C836.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C27CBF93-7E61-E111-AB30-003048F024C2.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/B46B6A67-8D61-E111-82C0-003048F118AA.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/36B93F53-7F61-E111-AA27-002481E0D524.root'
    ),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )




process.load('EGamma.EGammaAnalysisTools.electronIsoProducer_cfi')
process.isoPf = cms.Path(  process.elePFIso )

process.elePFIso.verbose = True


#my analyzer
#process.demo = cms.EDAnalyzer("ElectronAnalyzer")

#process.TFileService = cms.Service(
#    "TFileService",
#    fileName = cms.string("myhisto.root")
#    )

#process.pAna = cms.Path(process.demo)

process.schedule = cms.Schedule(process.isoPf)





