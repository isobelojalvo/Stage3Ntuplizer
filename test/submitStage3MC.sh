#!/bin/sh                                                                                                                                                                                                                                                                                                                     
#voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TStage3RatesEffMC2.py > SUBStage3.py
cat submitSecondaryFiles.py >> SUBStage3.py

mkdir /data/ojalvo/ggHttMC-v1-SUBStage3/
#submit dir = /data/ojalvo/ggHttMC-SUBStage3/submit
#make dag dir
mkdir -p /data/ojalvo/ggHttMC-v1-SUBStage3/dags
mkdir -p /data/ojalvo/ggHttMC-v1-SUBStage3/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v1/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-v1-SUBStage3/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ggHttMC-v1-SUBStage3/submit --output-dag-file=/data/ojalvo/ggHttMC-v1-SUBStage3/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v1/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-v1-SUBStage3/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage3.py 
farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=ggHtt-76X-TSG-MINI.txt --submit-dir=/data/ojalvo/ggHttMC-v1-SUBStage3/submit --output-dag-file=/data/ojalvo/ggHttMC-v1-SUBStage3/dags/dag  ggHttMC-v1   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBStage3.py farmout=True
rm SUBStage3.py

