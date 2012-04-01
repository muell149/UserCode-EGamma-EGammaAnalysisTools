import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START42_V12::All'

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
#    '/store/relval/CMSSW_4_4_0_pre8/RelValSingleGammaPt35/GEN-SIM-RECO/START44_V3-v1/0030/76A5E848-72C7-E011-AEA1-002354EF3BE1.root'
    'file:/data/benedet/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_0.root'
  ),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )


# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#my analyzer
process.demo = cms.EDAnalyzer("ElectronAnalyzer")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myhisto.root")
    )

process.pAna = cms.Path(process.demo)

process.schedule = cms.Schedule(process.pAna)





