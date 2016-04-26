import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TRatesAndEff",
                                  rctSource               = cms.InputTag("simRctDigis"),
                                  #gctTauJetsSource        = cms.InputTag("simCaloStage1LegacyFormatDigis","tauJets"),
                                  #gctIsoTauJetsSource     = cms.InputTag("simCaloStage1LegacyFormatDigis","isoTauJets"),
                                  l1ExtraTauSource        = cms.InputTag("uct2016EmulatorDigis","Tau"),
                                  l1ExtraIsoTauSource     = cms.InputTag("uct2016EmulatorDigis","IsoTau"),
                                  recoTau                 = cms.InputTag("slimmedTaus"),
                                  #remove all possible muons
                                  recoTauDiscriminatorIso = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                                  recoTauDiscriminatorMu  = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
                                  vertices                = cms.InputTag("offlinePrimaryVertices"),
                                  folderName              = cms.untracked.string("firstFolder"),
                                  recoPtCut               = cms.double(4),
                                  pfCands                 = cms.InputTag("packedPFCandidates"),
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis')
)
