#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TStage3RatesEffMC.py > SUBStage3.py
cat submitSecondaryFiles.py >> SUBStage3.py

rm -r /data/ojalvo/QCD_Pt-15to3000MC-$1-SUBStage3/
mkdir /data/ojalvo/QCD_Pt-15to3000MC-$1-SUBStage3/
#submit dir = /data/ojalvo/ QCD_Pt-15to3000MC-SUBStage3/submit
#make dag dir
mkdir -p /data/ojalvo/QCD_Pt-15to3000MC-$1-SUBStage3/dags
mkdir -p /data/ojalvo/QCD_Pt-15to3000MC-$1-SUBStage3/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/ojalvo/ QCD_Pt-15to3000MC-$1-SUBStage3/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ QCD_Pt-15to3000MC-$1-SUBStage3/submit --output-dag-file=/data/ojalvo/ QCD_Pt-15to3000MC-$1-SUBStage3/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/ojalvo/ QCD_Pt-15to3000MC-$1-SUBStage3/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage3.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=QCD_Pt-15to3000-MINIAOD.txt \
--submit-dir=/data/ojalvo/QCD_Pt-15to3000MC-$1-SUBStage3/submit \
--output-dag-file=/data/ojalvo/QCD_Pt-15to3000MC-$1-SUBStage3/dags/dag \
QCD_Pt-15to3000MC-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage3.py farmout=True rates=False data=False jets=True $2

rm SUBStage3.py

