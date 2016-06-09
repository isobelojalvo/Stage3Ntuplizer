# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --era=Run2_2016 --geometry=Extended2016,Extended2016Reco --customise=L1Trigger/Configuration/customiseReEmul.L1TEventSetupForHF1x1TPs --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU --customise=L1Trigger/Configuration/customiseUtils.L1TTurnOffUnpackStage2GtGmtAndCalo --conditions=auto:run2_data -n 100 --data --no_exec --no_output --filein=/store/data/Run2015D/ZeroBias1/RAW/v1/000/256/843/00000/FE8AD1BB-D05E-E511-B3A7-02163E01276B.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAW2DIGI',eras.Run2_2016)

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.register('farmout',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'options to set up cfi it is able to submit to condor')
options.register('data',True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'option to switch between data and mc')
options.register('jets',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'option to switch between jets and taus')
options.register('rates',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'option to switch between rates and eff (no secondary file map)')
options.parseArguments()

if options.rates is False :
    if options.jets is False :
        from L1Trigger.Stage3Ntuplizer.ggH_TSG_SecondaryFiles_cfi import *
    else :
        from L1Trigger.Stage3Ntuplizer.QCD_PT_15to3000_cfi import *

def formatLumis(lumistring, run) :
    lumis = (lrange.split('-') for lrange in lumistring.split(','))
    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
    return ['-'.join(l) for l in runlumis]

if options.farmout is False :
    print 'Getting files for run %d...' % options.runNumber
    if len(options.inputFiles) is 0 and options.inputFileList is '' :
        inputFiles = util.getFilesForRun(options.runNumber, options.dataStream)
    elif len(options.inputFileList) > 0 :
        with open(options.inputFileList) as f :
            inputFiles = list((line.strip() for line in f))
    else :
        inputFiles = cms.untracked.vstring(options.inputFiles)
        if len(inputFiles) is 0 :
            raise Exception('No files found for dataset %s run %d' % (options.dataStream, options.runNumber))

print 'Ok, time to analyze'

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

from Configuration.AlCa.GlobalTag import GlobalTag
if options.data is True :
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
else :
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_v6', '')


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

from L1Trigger.Configuration.customiseReEmul import L1TEventSetupForHF1x1TPs,L1TReEmulFromRAW 
process = L1TEventSetupForHF1x1TPs(process)

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

if options.jets is True :
    process.load("L1Trigger.Stage3Ntuplizer.l1RatesEffStage3Jets_cfi")
else :
    process.load("L1Trigger.Stage3Ntuplizer.l1RatesEffStage2_cfi")

process.l1NtupleProducer.ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
process.l1NtupleProducer.hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.l1NtupleProducer.folderName = cms.untracked.string("Stage2Taus")

process.simDtTriggerPrimitiveDigis.digiTag = 'muonDTDigis'  
process.simCscTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag( 'muonCSCDigis', 
                                                                              'MuonCSCComparatorDigi')
process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer       = cms.InputTag( 'muonCSCDigis', 
                                                                              'MuonCSCWireDigi' )  
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

if options.farmout is True :
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

##customizations for getting or not getting secondary file map
if options.farmout is False and options.rates is False:
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(inputFiles),
                                secondaryFileNames = cms.untracked.vstring(secondaryMap[inputFiles[0]])
                                )

if options.farmout is False and options.rates is True:
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(inputFiles)
                                )

#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple.root")
)

process.p = cms.Path(process.RawToDigi*process.simHcalTriggerPrimitiveDigis * process.simCaloStage2Layer1Digis * process.simCaloStage2Digis * process.l1NtupleProducer)


#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("l1TFullEvent.root"),
#    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
#    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
#)

#process.e = cms.EndPath(process.out)
