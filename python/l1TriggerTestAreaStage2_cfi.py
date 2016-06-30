import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("triggerTestArea",
                                  rctSource               = cms.InputTag("simRctDigis"),
                                  l1ExtraTauSource        = cms.InputTag("uct2016EmulatorDigis","Tau"),
                                  l1ExtraIsoTauSource     = cms.InputTag("uct2016EmulatorDigis","IsoTau"),
                                  stage2TauSource        = cms.InputTag("simCaloStage2Digis"),
                                  recoTau                 = cms.InputTag("slimmedTaus"),
                                  recoJets                 = cms.InputTag("slimmedJets"),
                                  #remove all possible muons
                                  recoTauDiscriminatorIso = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                                  recoTauDiscriminatorMu  = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
                                  vertices                = cms.InputTag("offlinePrimaryVertices"),
                                  recoPtCut               = cms.double(4),
                                  pfCands                 = cms.InputTag("packedPFCandidates"),
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis'),
                                  folderName              = cms.untracked.string("Stage3Taus"),##not used
)
