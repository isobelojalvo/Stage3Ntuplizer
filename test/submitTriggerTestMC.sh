#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TriggerTestArea.py > SUB.py
cat submitSecondaryFiles.py >> SUB.py

rm -r /data/ojalvo/ggHttMC-taus-CT-v2-SUB/
mkdir /data/ojalvo/ggHttMC-taus-CT-v2-SUB/
#submit dir = /data/ojalvo/ggHttMC-taus-CT-v2-SUB/submit
#make dag dir
mkdir -p /data/ojalvo/ggHttMC-taus-CT-v2-SUB/dags
mkdir -p /data/ojalvo/ggHttMC-taus-CT-v2-SUB/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-taus-CT-v2-SUB/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ggHttMC-taus-CT-v2-SUB/submit --output-dag-file=/data/ojalvo/ggHttMC-taus-CT-v2-SUB/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-taus-CT-v2-SUB/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUB.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=ggHtt-76X-TSG-MINI.txt \
--submit-dir=/data/ojalvo/ggHttMC-taus-CT-v2-SUB/submit \
--output-dag-file=/data/ojalvo/ggHttMC-taus-CT-v2-SUB/dags/dag \
ggHttMC-taus-CT-v2  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUB.py farmout=True rates=False data=False

rm SUB.py

#rates=True data=False