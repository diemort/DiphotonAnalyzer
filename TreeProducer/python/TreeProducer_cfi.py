import FWCore.ParameterSet.Config as cms

treeProducer = cms.EDAnalyzer('TreeProducer',

    # general parameters
    isData = cms.bool(True),
    sqrtS = cms.double(13.e3),
    outputFilename = cms.string('output.root'),
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # input collections
    hltMenuLabel = cms.string('HLT'),
    metLabel = cms.InputTag('flashggMetsCorr'),
    diphotonLabel = cms.InputTag('flashggDiPhotons'),
    vertexLabel = cms.InputTag('offlineSlimmedPrimaryVertices'),
    #electronLabel = cms.InputTag('flashggElectrons'),
    electronLabel = cms.InputTag('flashggSelectedElectrons'),
    #muonLabel = cms.InputTag('flashggMuons'),
    muonLabel = cms.InputTag('flashggSelectedMuons'),
    jetLabel = cms.InputTag('flashggFinalJets'),
    beamSpotLabel = cms.InputTag('offlineBeamSpot'),

    # "tight" single/double photon selection
    minPtSinglePhoton = cms.double(50.),
    minR9SinglePhoton = cms.double(0.94),
    maxEtaSinglePhoton = cms.double(2.5),
    minMassDiPhoton = cms.double(500.),

    # totem RP information extraction
    totemRPTracksLabel = cms.InputTag('ctppsLocalTrackLiteProducer'),
    useXiInterpolation = cms.bool(True),
    xiInterpolationFile = cms.FileInPath('DiphotonAnalyzer/TreeProducer/data/ctpps_optics_9mar2017.root'),
    fillNumLUTFile = cms.FileInPath('DiphotonAnalyzer/TreeProducer/data/fill_run_lut_v2.dat'),
    alignmentLUTFile = cms.FileInPath('DiphotonAnalyzer/TreeProducer/data/alignment_collection_v2.out'),

    # generator-level input collections
    genPartLabel = cms.InputTag('flashggPrunedGenParticles'),
    genPhoLabel = cms.InputTag('flashggGenPhotons'),
    maxGenLevelDeltaR = cms.double(5.0),
)
