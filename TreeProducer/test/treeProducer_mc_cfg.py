import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/user/j/juwillia/public/forLaurent/AC_microAOD/microAOD_zeta1_e-13_zeta2_e-13.root',
    )
)

# Trigger
#from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltHighLevel.HLTPaths = ['HLT_DoublePhoton60*', 'HLT_DoublePhoton85*']
process.hltHighLevel.HLTPaths = ['HLT_DoublePhoton60*']
process.hltHighLevel.throw = cms.bool(False)

process.load('DiphotonAnalyzer.TreeProducer.treeProducer_cfi')

# set some parameters to the run
process.treeProducer.isData = cms.bool(False)
process.treeProducer.minPtSinglePhoton = cms.double(50.)
process.treeProducer.minMassDiPhoton = cms.double(350.)
process.treeProducer.minR9SinglePhoton = cms.double(0.)
process.treeProducer.triggersList = process.hltHighLevel.HLTPaths
#FIXME for Incl GGprocess.treeProducer.metLabel = cms.InputTag('flashggMets')
process.treeProducer.pileupMCFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_mc_2017_25ns_WinterMC.root')
#process.treeProducer.pileupMCFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_mc.root')
#process.treeProducer.pileupDataFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_data16BCG_PPSruns_v2.root')
process.treeProducer.pileupDataFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_data16BCG_PPSruns_99bins.root')

process.p = cms.Path(
    process.hltHighLevel*
    process.treeProducer
)
