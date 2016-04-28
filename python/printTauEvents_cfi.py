import FWCore.ParameterSet.Config as cms


printTauEvents = cms.EDAnalyzer("printTauEvents",
                                  recoTau                 = cms.InputTag("slimmedTaus"),
                                  #remove all possible muons
                                  recoTauDiscriminatorIso = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                                  recoTauDiscriminatorMu  = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
                                  vertices                = cms.InputTag("offlinePrimaryVertices"),
                                  recoPtCut               = cms.double(4),
                                  maxRecoPt               = cms.double(25)
#                                  pfCands                 = cms.InputTag("packedPFCandidates"),
)
