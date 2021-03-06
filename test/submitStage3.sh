#!/bin/sh                                                                                                                              
#voms-proxy-init --voms cms --valid 100:00                                                                                              

cat testL1TStage3RatesEffMC.py > SUBStage3.py
cat submit.py >> SUBStage3.py

rm -r /data/ojalvo/ZeroBias-259721-$1-SUBStage3/
mkdir /data/ojalvo/ZeroBias-259721-$1-SUBStage3/
#submit dir = /data/ojalvo/ZeroBias-259721-$1-SUBStage3/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-259721-$1-SUBStage3/dags
mkdir -p /data/ojalvo/ZeroBias-259721-$1-SUBStage3/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-259721-$1-SUBStage3/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-259721-$1-SUBStage3/submit --output-dag-file=/data/ojalvo/ZeroBias-259721-$1-SUBStage3/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-259721-$1-SUBStage3/ ZeroBias-259721-$1-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBStage3.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-259721-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-259721-$1-SUBStage3/submit \
--output-dag-file=/data/ojalvo/ZeroBias-259721-$1-SUBStage3/dags/dag \
ZeroBias-259721-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage3.py farmout=True rates=True data=True



rm -r /data/ojalvo/ZeroBias-260627-$1-SUBStage3/
mkdir /data/ojalvo/ZeroBias-260627-$1-SUBStage3/
#submit dir = /data/ojalvo/ZeroBias-260627-$1-SUBStage3/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-260627-$1-SUBStage3/dags
mkdir -p /data/ojalvo/ZeroBias-260627-$1-SUBStage3/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-$1-SUBStage3/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-260627-$1-SUBStage3/submit --output-dag-file=/data/ojalvo/ZeroBias-260627-$1-SUBStage3/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-$1-SUBStage3/ ZeroBias-260627-$1   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBStage3.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-260627-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-260627-$1-SUBStage3/submit \
--output-dag-file=/data/ojalvo/ZeroBias-260627-$1-SUBStage3/dags/dag \
ZeroBias-260627-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage3.py farmout=True rates=True data=True

rm SUBStage3.py

