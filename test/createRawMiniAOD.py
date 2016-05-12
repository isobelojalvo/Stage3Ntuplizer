import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

import EventFilter.L1TXRawToDigi.util as util

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputSecondaryFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('outputFileName', 'file', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual outputFile Name')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.parseArguments()

def formatLumis(lumistring, run) :
    lumis = (lrange.split('-') for lrange in lumistring.split(','))
    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
    return ['-'.join(l) for l in runlumis]

print 'Getting files for run %d...' % options.runNumber
if len(options.inputFiles) is 0 and options.inputFileList is '' :
    inputFiles = util.getFilesForRun(options.runNumber, options.dataStream)
elif len(options.inputFileList) > 0 :
    with open(options.inputFileList) as f :
        inputFiles = list((line.strip() for line in f))
else :
    inputFiles = cms.untracked.vstring(options.inputFiles)

if len(options.inputSecondaryFiles) is 0 and options.inputFileList is '' :
    inputSecondaryFiles = util.getFilesForRun(options.runNumber, options.dataStream)
elif len(options.inputFileList) > 0 :
    with open(options.inputFileList) as f :
        inputSecondaryFiles = list((line.strip() for line in f))
else :
    inputSecondaryFiles = cms.untracked.vstring(options.inputSecondaryFiles)
if len(inputFiles) is 0 :
    raise Exception('No files found for dataset %s run %d' % (options.dataStream, options.runNumber))
print 'Ok, time to analyze'

secondaryMap = {
    "/store/mc/RunIIFall15DR76/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/003220F0-CDA8-E511-817B-5065F3816291.root": ["/store/mc/RunIIFall15MiniAODv1/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/128A0602-93A9-E511-9262-A0000420FE80.root","/store/mc/RunIIFall15MiniAODv2/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/92A604E8-CFB9-E511-AF99-A0000420FE80.root"],
    "root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/MINIAOD/16Dec2015-v1/10000/0A9C7582-B8A8-E511-94C9-0CC47A78A478.root": [
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/08D44671-7082-E511-8851-02163E0141EA.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/288AA947-7082-E511-8314-02163E013909.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/4606E000-7382-E511-B468-02163E01476B.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/586161DC-7182-E511-A51B-02163E014281.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/5C067250-7182-E511-9AE9-02163E011B44.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/7CBE3F06-6182-E511-A818-02163E0138E1.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/C05D9969-7082-E511-A1CF-02163E014529.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/C0A875B0-6182-E511-BFAB-02163E01464C.root'],
    "file:/data/ojalvo/localFiles/Run2015D/SingleMuon/MINIAOD/16Dec2015-v1/10000/0A9C7582-B8A8-E511-94C9-0CC47A78A478.root":[
        'file:/data/ojalvo/localFiles/Run2015D/SingleMuon/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/4606E000-7382-E511-B468-02163E01476B.root'],
    "file:l1TFullEvent-Mini-RECO-260267-LUMI-1643.root":[
        "file:/data/ojalvo/localFiles/Run2015D/SingleMuon/RAW/v1/000/260/627/00000/7CBE3F06-6182-E511-A818-02163E0138E1.root"],
    "file:pickevents-highPtTaus2.root":[
        'root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/0A6BACC1-A79C-E511-8606-0CC47A4C8EA8.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/129E931E-B09C-E511-93F4-0CC47A4D76BE.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/2210B287-AB9C-E511-A185-0CC47A4C8F08.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/4A51B3AC-A89C-E511-9D46-0CC47A4C8EA8.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/4A8DA380-A89C-E511-B535-0CC47A4D767E.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/4CD53D8A-AB9C-E511-BFD7-0026189438E8.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/526D5BBE-B19C-E511-944F-0025905A60CA.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/76E9CDE3-B09C-E511-89FE-0CC47A4D76AA.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/86E0E6B3-B19C-E511-A614-0025905A6070.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/9A2BE4C4-B29C-E511-822B-0025905A6094.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/A686A218-AD9C-E511-A313-0CC47A4C8EA8.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/AE68D4BE-B39C-E511-9C34-0025905A60CE.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/AEC19F8E-B59C-E511-922E-0CC47A4D764C.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/B0C09325-AC9C-E511-A90F-0025905B8568.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/B0C650DB-A69C-E511-9331-0CC47A4C8EA8.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/B266BEBF-B29C-E511-9187-0CC47A4D76BE.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/DE6A120B-A29C-E511-BBE8-00261894397F.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/ECC3251C-AD9C-E511-871C-002590593878.root','root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v12-v1/50000/FE388856-AC9C-E511-A0F5-0CC47A4C8E0E.root'],
    "file:pickEventsHighL1LowReco-MINI.root":[
        "file:pickEventsHighL1LowReco.root"
        ],
    "file:pickEvent-769674-MINI.root":[
        "file:pickEvent-769674-RAW.root"
        ],
    "file:pickEvent-547237-MINI.root":[
        "file:pickEvent-547237-RAW.root"
        ]
}

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
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(inputFiles),
                            #secondaryFileNames = cms.untracked.vstring(secondaryMap[options.inputFiles[0]]),
                            secondaryFileNames = cms.untracked.vstring(inputSecondaryFiles),
                            #secondaryFileNames = cms.untracked.vstring("file:/data/ojalvo/localFiles/Run2015D/SingleMuon/RECO/16Dec2015-v1/10016/FAF49A05-DCA7-E511-8A3A-0CC47A4D7614.root")
)

outputFile = '/data/' + os.environ['USER'] + '/l1tCaloSummary-' + str(options.runNumber) + '.root'

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFileName),
    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)


#Output
#process.TFileService = cms.Service(
#	"TFileService",
#	fileName = cms.string("event-ggH-546248-Full.root")
#)

#process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

process.e = cms.EndPath(process.out)
