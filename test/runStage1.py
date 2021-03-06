# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: STAGE1 --python_filename=runStage1.py -s RAW2DIGI --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --datatier FEVT --era=Run2_25ns --conditions=auto:run2_data --data --no_exec -n 1000 --filein=/store/data/Run2015D/JetHT/RAW/v1/000/260/627/00000/0021ABE7-3282-E511-B229-02163E012312.root
import os
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('RAW2DIGI',eras.Run2_25ns)

import EventFilter.L1TXRawToDigi.util as util

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
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
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
    annotation = cms.untracked.string('STAGE1 nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

#process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('FEVT'),
#        filterName = cms.untracked.string('')
#    ),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    fileName = cms.untracked.string('STAGE1_RAW2DIGI.root'),
#    outputCommands = process.RECOSIMEventContent.outputCommands,
#    splitLevel = cms.untracked.int32(0)
#)

# Additional output definition


# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step,process.RECOSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW 

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
    cms.InputTag('hcalDigis'),
    cms.InputTag('hcalDigis')
    )
#process.L1TReEmul = cms.Sequence(process.simHcalTriggerPrimitiveDigis * process.SimL1Emulator)
process.simDtTriggerPrimitiveDigis.digiTag = 'muonDTDigis'  
process.simCscTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag( 'muonCSCDigis', 'MuonCSCComparatorDigi')
process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer       = cms.InputTag( 'muonCSCDigis', 'MuonCSCWireDigi' )  


####
process.simRctDigis.ecalDigis = cms.VInputTag( cms.InputTag( 'ecalDigis:EcalTriggerPrimitives' ) )
process.simRctDigis.hcalDigis = cms.VInputTag('simHcalTriggerPrimitiveDigis')
process.simRpcTriggerDigis.label         = 'muonRPCDigis'
process.simRpcTechTrigDigis.RPCDigiLabel  = 'muonRPCDigis'

if options.jets is True :
    process.load("L1Trigger.Stage3Ntuplizer.l1RatesEffStage2Jets_cfi")
    process.l1NtupleProducer.folderName = cms.untracked.string("Stage1Jets")
    process.l1NtupleProducer.stage2JetSource = cms.InputTag("garbage")
    process.l1NtupleProducer.stage1JetSource    = cms.InputTag("simCaloStage1FinalDigis","preGtJets")
else :
    process.load("L1Trigger.Stage3Ntuplizer.l1RatesEffStage2_cfi")
    process.l1NtupleProducer.ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
    process.l1NtupleProducer.hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
    process.l1NtupleProducer.folderName = cms.untracked.string("Stage1Taus")
#process.l1NtupleProducer.stage2TauSource = cms.InputTag("simCaloStage1FinalDigis","rlxTaus")

#process.L1TReEmulPath = cms.Path(process.L1TReEmul)

#process = L1TReEmulFromRAW(process)

# End of customisation functions

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
    fileName = cms.string("l1TNtuple2.root")
    )

#if options.data is True :
#    process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)
#else :
#    process.p = cms.Path(process.RawToDigi*process.L1TReEmul*process.uct2016EmulatorDigis*process.l1NtupleProducer)
process.p = cms.Path(process.RawToDigi * process.simHcalTriggerPrimitiveDigis * process.SimL1Emulator * process.l1NtupleProducer)


#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("l1TFullEvent.root"),
#    outputCommands = cms.untracked.vstring('keep *_*_*_RAW2DIGI') #'keep *_*_*_L1TCaloSummaryTest')
##    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
#)
#process.e = cms.EndPath(process.out)
