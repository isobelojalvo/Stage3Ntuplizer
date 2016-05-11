#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat runStage1.py > SUBStage1.py
cat submitSecondaryFiles.py >> SUBStage1.py

rm -r /data/ojalvo/ggHttMC-SUBStage1/
mkdir /data/ojalvo/ggHttMC-SUBStage1/
#submit dir = /data/ojalvo/ggHttMC-SUBStage1/submit
#make dag dir
mkdir -p /data/ojalvo/ggHttMC-SUBStage1/dags
mkdir -p /data/ojalvo/ggHttMC-SUBStage1/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-SUBStage1/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ggHttMC-SUBStage1/submit --output-dag-file=/data/ojalvo/ggHttMC-SUBStage1/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-SUBStage1/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage1.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=ggHtt-76X-TSG-MINI.txt \
--submit-dir=/data/ojalvo/ggHttMC-SUBStage1/submit \
--output-dag-file=/data/ojalvo/ggHttMC-SUBStage1/dags/dag \
ggHttMC  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage1.py farmout=True rates=False data=False

rm SUBStage1.py

