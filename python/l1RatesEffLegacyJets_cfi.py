import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TRatesAndEffJets",
                                  l1ExtraJetSource   = cms.InputTag("uct2016EmulatorDigis","Central"),
                                  stage2JetSource    = cms.InputTag("caloStage2Digis","Jet"),
#                                  stage2TauSource    = cms.InputTag("caloStage2Digis","Tau"),
                                  stage1JetSource    = cms.InputTag("simCaloStage1FinalDigis","preGtJets"),
                                  recoJets           = cms.InputTag("slimmedJets"),
                                  recoJetsAK8        = cms.InputTag("slimmedJetsAK8"),
                                  vertices           = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  recoPtCut          = cms.double(15),
                                  pfCands            = cms.InputTag("packedPFCandidates"),
                                  ecalDigis          = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis          = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  folderName         = cms.untracked.string("Stage3Jets"),
)
