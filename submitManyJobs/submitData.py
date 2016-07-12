process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        $inputFileNames
                                      )
    )
import PhysicsTools.PythonAnalysis.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt').getVLuminosityBlockRange()

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("$outputFileName")
)
