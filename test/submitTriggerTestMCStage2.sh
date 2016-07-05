#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TriggerTestAreaStage2.py > SUB-Stage2.py
cat submitSecondaryFiles.py >> SUB-Stage2.py

rm -r /data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/
mkdir /data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/
#submit dir = /data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/submit
#make dag dir
mkdir -p /data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/dags
mkdir -p /data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/submit --output-dag-file=/data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUB-Stage2.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=QCD_Pt-15to3000-MINIAOD.txt \
--submit-dir=/data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/submit \
--output-dag-file=/data/ojalvo/QCD_Pt-15to3000-triggerTestArea-taus-SUB-Stage2/dags/dag \
QCD_Pt-15to3000-taus  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUB-Stage2.py farmout=True rates=False data=False

rm SUB-Stage2.py

#rates=True data=False