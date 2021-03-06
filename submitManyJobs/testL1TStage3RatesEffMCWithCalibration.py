import os
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1TCaloSummaryTest")

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
options.register('rates',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'option to switch between rates and eff (no secondary file map)')
options.register('jets',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'option to switch between jets and taus')
options.register('tauIsolationFactor', 0.3, VarParsing.multiplicity.singleton, VarParsing.varType.float, 'tau relative isolation factor')
options.register('tauMiscActivityThresh', 0.25, VarParsing.multiplicity.singleton, VarParsing.varType.float, 'tau relative isolation factor')#options.tauMiscActivityThresh
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
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
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


# To get L1 CaloParams
process.load('L1Trigger.L1TCalorimeter.caloStage2Params_2016_v2_1_cfi')
# To get CaloTPGTranscoder
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)


from L1Trigger.Configuration.customiseReEmul import L1TEventSetupForHF1x1TPs,L1TReEmulFromRAW 
process = L1TEventSetupForHF1x1TPs(process)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
    cms.InputTag('hcalDigis'),
    cms.InputTag('hcalDigis')
    )
process.L1TReEmul = cms.Sequence(process.simHcalTriggerPrimitiveDigis * process.SimL1Emulator)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

if options.jets is True :
    process.load("L1Trigger.Stage3Ntuplizer.l1RatesEffStage3Jets_cfi")
else :
    process.load("L1Trigger.Stage3Ntuplizer.l1RatesEffStage3_cfi")
    process.l1NtupleProducer.ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
    process.l1NtupleProducer.hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
    process.l1NtupleProducer.stage2TauSource = cms.InputTag("garbage")
    process.l1NtupleProducer.stage1TauSource = cms.InputTag("garbage")
    process.l1NtupleProducer.stage1IsoTauSource = cms.InputTag("garbage")

process.uct2016EmulatorDigis.useECALLUT = cms.bool(True)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(True)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.tauSeed = cms.uint32(10)
process.uct2016EmulatorDigis.tauIsolationFactor = cms.double(options.tauIsolationFactor)
process.uct2016EmulatorDigis.tauMiscActivityThreshold = cms.double(options.tauMiscActivityThresh)
###
process.uct2016EmulatorDigis.tauLUT00 = cms.vdouble(1.7)#0-10GeV        
process.uct2016EmulatorDigis.tauLUT01 = cms.vdouble(1.7)#10-16GeV       
process.uct2016EmulatorDigis.tauLUT02 = cms.vdouble(1.5)#16-18GeV       
process.uct2016EmulatorDigis.tauLUT03 = cms.vdouble(1.3)#18-20GeV       
process.uct2016EmulatorDigis.tauLUT04 = cms.vdouble(1.15)#20-22GeV      
process.uct2016EmulatorDigis.tauLUT05 = cms.vdouble(1.12)#22-25GeV      
process.uct2016EmulatorDigis.tauLUT06 = cms.vdouble(1.08)#25-28GeV      
process.uct2016EmulatorDigis.tauLUT07 = cms.vdouble(1.07)#28-30GeV      
process.uct2016EmulatorDigis.tauLUT08 = cms.vdouble(1.06)#30-32GeV      
process.uct2016EmulatorDigis.tauLUT09 = cms.vdouble(1.05)#32-34GeV      
process.uct2016EmulatorDigis.tauLUT10 = cms.vdouble(1.02)#34-40GeV      
process.uct2016EmulatorDigis.tauLUT11 = cms.vdouble(1.02)#40-50GeV      
process.uct2016EmulatorDigis.tauLUT12 = cms.vdouble(1.01)#50-60GeV      
process.uct2016EmulatorDigis.tauLUT13 = cms.vdouble(1.01)#60-70GeV      
process.uct2016EmulatorDigis.tauLUT14 = cms.vdouble(1.01)#70-80GeV      
process.uct2016EmulatorDigis.tauLUT15 = cms.vdouble(1.01)#80-90GeV      
process.uct2016EmulatorDigis.tauLUT16 = cms.vdouble(1.01)#90-100GeV     
process.uct2016EmulatorDigis.tauLUT17 = cms.vdouble(1.01)#max is 125 GeV

###
print "tau Isolation Factor: " 
print options.tauIsolationFactor

#unpack the data or use the readout
if options.data is True :
    process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
    process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")
else :
    process.l1NtupleProducer.ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
    process.l1NtupleProducer.hcalDigis = cms.InputTag("hcalDigis")
    process.uct2016EmulatorDigis.ecalToken = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
    process.uct2016EmulatorDigis.hcalToken = cms.InputTag("hcalDigis")
    
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
print "now running through"
if options.data is True :
    process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)
else :
    process.p = cms.Path(process.RawToDigi*process.L1TReEmul*process.uct2016EmulatorDigis*process.l1NtupleProducer)

#process.p = cms.Path(process.L1TReEmul*process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())

#process.e = cms.EndPath(process.out)
