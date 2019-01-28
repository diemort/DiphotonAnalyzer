import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/data/Run2016B/MinimumBias/AOD/23Sep2016-v1/50000/040C2700-A388-E611-AF56-00259073E390.root',
#'/store/data/Run2016C/L1MinimumBias/AOD/PromptReco-v2/000/275/656/00000/3AD13D68-3F3B-E611-8BC1-02163E0143EF.root'
#'/store/data/Run2016C/DoubleEG/AOD/07Aug17-v1/50000/000F79D0-929E-E711-8092-0CC47A4D76B6.root'
#'/store/data/Run2016C/DoubleEG/AOD/23Sep2016-v1/120000/6A69C51A-10B8-E611-88DB-0CC47A78A4A6.root'
#'/store/data/Run2016C/DoubleEG/AOD/23Sep2016-v1/50000/0089FDD7-1783-E611-A026-90B11C2AA16C.root'
        '/store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/00000/D0A00A92-EB8A-E811-899C-008CFA1979AC.root',
    )
)
# Trigger

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltHighLevel.HLTPaths = ['HLT_L1MinimumBias*']
process.hltHighLevel.HLTPaths = ['*']

process.load('DiphotonAnalyzer.TreeProducer.mbTreeProducer_cfi')

# set some parameters to the run

process.p = cms.Path(
    process.hltHighLevel*
    process.mbTreeProducer
)
