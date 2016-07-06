#!/bin/sh

voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TriggerTestArea.py > TT-SUB-Stage3.py
cat submitSecondaryFiles.py >> TT-SUB-Stage3.py

rm -r /data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/
mkdir /data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/
#submit dir = /data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/submit
#make dag dir
mkdir -p /data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/dags
mkdir -p /data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/submit --output-dag-file=/data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/TT-SUB-Stage3.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=QCD_Pt-15to3000-MINIAOD.txt \
--submit-dir=/data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/submit \
--output-dag-file=/data/ojalvo/QCD_Pt-15to3000-taus-TT-SUB-Stage3/dags/dag \
QCD_Pt-15to3000-taus  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/TT-SUB-Stage3.py farmout=True rates=False data=False jets=True

rm TT-SUB-Stage3.py

#rates=True data=False