#!/bin/sh                                                                                                                              
voms-proxy-init --voms cms --valid 100:00                                                                                              

cat runStage1.py > SUBStage1.py
cat submit.py >> SUBStage1.py

rm -r /data/ojalvo/ZeroBias-260627-SUBStage1/
mkdir /data/ojalvo/ZeroBias-260627-SUBStage1/
#submit dir = /data/ojalvo/ZeroBias-260627-SUBStage1/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-260627-SUBStage1/dags
mkdir -p /data/ojalvo/ZeroBias-260627-SUBStage1/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-SUBStage1/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-260627-SUBStage1/submit --output-dag-file=/data/ojalvo/ZeroBias-260627-SUBStage1/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-SUBStage1/ ZeroBias-260627-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBStage1.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-260627-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-260627-SUBStage1/submit \
--output-dag-file=/data/ojalvo/ZeroBias-260627-SUBStage1/dags/dag \
ZeroBias-260627  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage1.py farmout=True rates=True data=True

rm SUBStage1.py

