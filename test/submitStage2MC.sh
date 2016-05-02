#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TStage2RatesEff.py > SUBStage2.py
cat submitSecondaryFiles.py >> SUBStage2.py

rm -r /data/ojalvo/ggHttMC-SUBStage2/
mkdir /data/ojalvo/ggHttMC-SUBStage2/
#submit dir = /data/ojalvo/ggHttMC-SUBStage2/submit
#make dag dir
mkdir -p /data/ojalvo/ggHttMC-SUBStage2/dags
mkdir -p /data/ojalvo/ggHttMC-SUBStage2/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-SUBStage2/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ggHttMC-SUBStage2/submit --output-dag-file=/data/ojalvo/ggHttMC-SUBStage2/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-SUBStage2/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage2.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=ggHtt-76X-TSG-MINI.txt \
--submit-dir=/data/ojalvo/ggHttMC-SUBStage2/submit \
--output-dag-file=/data/ojalvo/ggHttMC-SUBStage2/dags/dag \
ggHttMC  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage2.py farmout=True rates=False data=False

rm SUBStage2.py

