import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")

import EventFilter.L1TXRawToDigi.util as util

from FWCore.ParameterSet.VarParsing import VarParsing
from L1Trigger.Stage3Ntuplizer.ggH_TSG_SecondaryFiles_cfi import *
#from L1Trigger.Stage3Ntuplizer.HTT_AODSIM_cfi import *

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.register('farmout',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'options to set up cfi it is able to submit to condor')
options.register('data',True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'option to switch between data and mc')
options.parseArguments()

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

#secondaryMap = {
#    "root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/10000/04EDCDA3-916F-E511-AD11-0025905938B4.root": [
#        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/257/487/00000/28774BB9-EF66-E511-A328-02163E011CC3.root',
#        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/257/487/00000/8A42A2DA-EF66-E511-AC62-02163E012988.root',
#        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/257/487/00000/C24B52CB-EF66-E511-969B-02163E0119F6.root']
#}

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
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_ForTSGStudies_fromRunIIFall15DR76_v1', '')

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Stage3Ntuplizer.printTauEvents_cfi")
process.printTauEvents.recoPtCut = cms.double(35)
process.printTauEvents.maxRecoPt = cms.double(100)

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

if options.farmout is False :
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(400000) )
else :
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
 
if options.farmout is False :
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(inputFiles),
                                #secondaryFileNames = cms.untracked.vstring(secondaryMap[options.inputFiles[0]])
                                )
    #Output
    process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple2.root")
        )



process.p = cms.Path(process.printTauEvents)

#process.e = cms.EndPath(process.out)
