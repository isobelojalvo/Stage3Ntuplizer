import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TRatesAndEff",
                                  rctSource               = cms.InputTag("simRctDigis"),
                                  l1ExtraTauSource        = cms.InputTag("uct2016EmulatorDigis","Tau"),
                                  l1ExtraIsoTauSource     = cms.InputTag("uct2016EmulatorDigis","IsoTau"),
                                  stage2TauSource        = cms.InputTag("simCaloStage2Digis"),
                                  stage1TauSource        = cms.InputTag("simCaloStage1FinalDigis","rlxTaus"),
                                  stage1IsoTauSource     = cms.InputTag("simCaloStage1FinalDigis","isoTaus"),
                                  recoTau                 = cms.InputTag("slimmedTaus"),
                                  #remove all possible muons
                                  recoTauDiscriminatorIso = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                                  recoTauDiscriminatorMu  = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
                                  vertices                = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  recoPtCut               = cms.double(4),
                                  pfCands                 = cms.InputTag("packedPFCandidates"),
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis'),
                                  folderName              = cms.untracked.string("Stage2Taus")
)
