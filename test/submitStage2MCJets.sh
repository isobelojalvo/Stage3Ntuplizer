#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TStage2RatesEff.py > SUBStage2.py
cat submitSecondaryFiles.py >> SUBStage2.py

rm -r /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/
mkdir /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/
#submit dir = /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/submit
#make dag dir
mkdir -p /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/dags
mkdir -p /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/submit --output-dag-file=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/ ZeroBias-560627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage2.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=QCD_Pt-15to3000-MINIAOD.txt \
--submit-dir=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/submit \
--output-dag-file=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBStage2/dags/dag \
QCD_Pt-15to3000MC-Jets-5  \
--vsize-limit=7000 \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage2.py farmout=True rates=False data=False jets=True

rm SUBStage2.py

