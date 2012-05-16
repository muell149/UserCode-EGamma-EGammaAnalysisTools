import FWCore.ParameterSet.Config as cms

process = cms.Process("ExISO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/relval/CMSSW_5_2_5/RelValSingleElectronPt35/GEN-SIM-RECO/PU_START52_V9_PU_E7TeV_FIX_2_BX156_N20_BX50_Bm8_2_special_120508-v2/0010/0021B7A8-5C99-E111-BD90-003048D37560.root'
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C61C71DC-8061-E111-AAEB-0025B32036D2.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C49AC4E5-7F61-E111-8FB2-003048F1C836.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C27CBF93-7E61-E111-AB30-003048F024C2.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/B46B6A67-8D61-E111-82C0-003048F118AA.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/36B93F53-7F61-E111-AA27-002481E0D524.root'
    )
    )




process.load('EGamma.EGammaAnalysisTools.electronIsoProducer_cfi')


process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myIsohisto.root")
    )
process.eleIsoAnalyzer = cms.EDAnalyzer("ElectronIsoAnalyzer",
                                        verbose = cms.untracked.bool(True),
                                        Electrons = cms.InputTag('gsfElectrons'),
                                        IsoValElectrons = cms.VInputTag(cms.InputTag('elePFIso:chIsoForGsfEle'),
                                                                        cms.InputTag('elePFIso:phIsoForGsfEle'),
                                                                        cms.InputTag('elePFIso:nhIsoForGsfEle'))
                                        )




process.p = cms.Path(  process.elePFIso+process.eleIsoAnalyzer )
process.elePFIso.verbose = True


