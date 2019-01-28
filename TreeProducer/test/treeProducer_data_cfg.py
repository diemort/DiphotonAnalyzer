import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/group/phys_pps/diphoton/DoubleEG/DoubleEG_microAOD_ReReco_run2016B/181201_222130/0001/myMicroAODOutputFile_1033.root',
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

# Trigger
#from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltHighLevel.HLTPaths = ['HLT_DoublePhoton60_*', 'HLT_DoublePhoton85_*']
process.hltHighLevel.throw = cms.bool(False)

# Proton reconstruction
process.load('RecoCTPPS.ProtonReconstruction.year_2016.ctppsProtonReconstruction_cfi')

# set some parameters to the run
process.load('DiphotonAnalyzer.TreeProducer.treeProducer_cfi')
process.treeProducer.minPtSinglePhoton = cms.double(50.)
process.treeProducer.minMassDiPhoton = cms.double(350.)
process.treeProducer.minR9SinglePhoton = cms.double(0.)
process.treeProducer.triggersList = process.hltHighLevel.HLTPaths

process.p = cms.Path(
    process.hltHighLevel*
    process.ctppsProtonReconstruction*
    process.treeProducer
)
