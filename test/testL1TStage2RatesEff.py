# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --era=Run2_2016 --geometry=Extended2016,Extended2016Reco --customise=L1Trigger/Configuration/customiseReEmul.L1TEventSetupForHF1x1TPs --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU --customise=L1Trigger/Configuration/customiseUtils.L1TTurnOffUnpackStage2GtGmtAndCalo --conditions=auto:run2_data -n 100 --data --no_exec --no_output --filein=/store/data/Run2015D/ZeroBias1/RAW/v1/000/256/843/00000/FE8AD1BB-D05E-E511-B3A7-02163E01276B.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAW2DIGI',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2015D/ZeroBias1/RAW/v1/000/256/843/00000/FE8AD1BB-D05E-E511-B3A7-02163E01276B.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TEventSetupForHF1x1TPs,L1TReEmulFromRAW 

#call to customisation function L1TEventSetupForHF1x1TPs imported from L1Trigger.Configuration.customiseReEmul
process = L1TEventSetupForHF1x1TPs(process)

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
#process = L1TReEmulFromRAW(process)
process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
    cms.InputTag('hcalDigis'),
    cms.InputTag('hcalDigis')
    )


from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import *
from L1Trigger.L1TCaloLayer1.simCaloStage2Layer1Digis_cfi import simCaloStage2Layer1Digis
from L1Trigger.L1TCalorimeter.simCaloStage2Digis_cfi import simCaloStage2Digis

process.simDtTriggerPrimitiveDigis.digiTag = 'muonDTDigis'  
process.simCscTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag( 'muonCSCDigis', 'MuonCSCComparatorDigi')
process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer       = cms.InputTag( 'muonCSCDigis', 'MuonCSCWireDigi' )  

process.simTwinMuxDigis.RPC_Source         = cms.InputTag('muonRPCDigis')
# When available, this will switch to TwinMux input Digis:
process.simTwinMuxDigis.DTDigi_Source      = cms.InputTag("dttfDigis")
process.simTwinMuxDigis.DTThetaDigi_Source = cms.InputTag("dttfDigis")
process.simOmtfDigis.srcRPC                = cms.InputTag('muonRPCDigis')
process.simBmtfDigis.DTDigi_Source         = cms.InputTag("simTwinMuxDigis")
process.simBmtfDigis.DTDigi_Theta_Source   = cms.InputTag("dttfDigis")
process.simEmtfDigis.CSCInput              = cms.InputTag("csctfDigis")
process.simOmtfDigis.srcCSC                = cms.InputTag("csctfDigis")
process.simCaloStage2Layer1Digis.ecalToken = cms.InputTag("ecalDigis","EcalTriggerPrimitives")

#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple.root")
)

process.p = cms.Path(process.RawToDigi*process.simHcalTriggerPrimitiveDigis * process.simCaloStage2Layer1Digis * process.simCaloStage2Digis )


#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("l1TFullEvent.root"),
#    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
#    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
#)

#process.e = cms.EndPath(process.out)
