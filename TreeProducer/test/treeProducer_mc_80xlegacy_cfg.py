import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
#'/store/group/phys_pps/diphoton/VBFHToGG_M125_13TeV_amcatnlo_pythia8/lforthom-microAOD-ctpps_VBFHToGG_M125_13TeV_amcatnlo_pythia8_v1/170315_142732/0000/myMicroAODOutputFile_61.root',
#'/store/group/phys_pps/diphoton/VBFHToGG_M125_13TeV_amcatnlo_pythia8/lforthom-microAOD-ctpps_VBFHToGG_M125_13TeV_amcatnlo_pythia8_v1/170315_142732/0000/myMicroAODOutputFile_62.root',
#'/store/group/phys_pps/diphoton/VBFHToGG_M125_13TeV_amcatnlo_pythia8/lforthom-microAOD-ctpps_VBFHToGG_M125_13TeV_amcatnlo_pythia8_v1/170315_142732/0000/myMicroAODOutputFile_63.root'
#'/store/group/phys_pps/diphoton/GammaGammaToEE_13TeV_lpair/myMicroAODOutputFile_GammaGammaEE_lpair_elastic.root'
#'/store/group/phys_pps/diphoton/GammaGammaToEE_13TeV_lpair/myMicroAODOutputFile_GammaGammaEE_lpair_singleinelastic.root'
#'/store/group/phys_pps/diphoton/GammaGammaToEE_13TeV_lpair/myMicroAODOutputFile_GammaGammaEE_lpair_doubleinelastic.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180606_162427/0000/myMicroAODOutputFile_99.root',
'/store/group/phys_higgs/cmshgg/lforthom/flashgg/RunIIPPS16_94X_reminiAOD/2_7_6/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/RunIIPPS16_94X_reminiAOD-2_7_6-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/181217_204547/0000/myMicroAODOutputFile_1.root',
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_1/3_1_1/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/RunIIFall17-3_1_1-3_1_1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/180828_081414/0000/myMicroAODOutputFile_610.root',
    )
)

# Trigger
#from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltHighLevel.HLTPaths = ['HLT_DoublePhoton60*', 'HLT_DoublePhoton85*']
process.hltHighLevel.throw = cms.bool(False)

process.load('DiphotonAnalyzer.TreeProducer.treeProducer_cfi')

# set some parameters to the run
process.treeProducer.isData = cms.bool(False)
process.treeProducer.minPtSinglePhoton = cms.double(50.)
process.treeProducer.minMassDiPhoton = cms.double(350.)
process.treeProducer.minR9SinglePhoton = cms.double(0.)
process.treeProducer.triggersList = process.hltHighLevel.HLTPaths
###process.treeProducer.metLabel = cms.InputTag('flashggMets')##FIXME
process.treeProducer.pileupMCFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_mc_2016_25ns_Moriond17MC.root')
#process.treeProducer.pileupMCFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_mc.root')
#process.treeProducer.pileupDataFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_data16BCG_PPSruns_v2.root')
process.treeProducer.pileupDataFile = cms.untracked.FileInPath('DiphotonAnalyzer/TreeProducer/data/pileup_data16BCG_PPSruns_75bins.root')

process.p = cms.Path(
    process.hltHighLevel*
    process.treeProducer
)
