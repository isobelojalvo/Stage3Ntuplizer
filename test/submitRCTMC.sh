#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   

cat testL1TRCTRatesEff.py > SUBRct.py
cat submitSecondaryFiles.py >> SUBRct.py

rm -r /data/ojalvo/ggHttMC-$1-SUBRct/
mkdir /data/ojalvo/ggHttMC-$1-SUBRct/
#submit dir = /data/ojalvo/ggHttMC-$1-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/ggHttMC-$1-SUBRct/dags
mkdir -p /data/ojalvo/ggHttMC-$1-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-$1-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ggHttMC-$1-SUBRct/submit --output-dag-file=/data/ojalvo/ggHttMC-$1-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ggHttMC-$1-SUBRct/ ZeroBias-260627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=ggHtt-76X-TSG-MINI.txt \
--submit-dir=/data/ojalvo/ggHttMC-$1-SUBRct/submit \
--output-dag-file=/data/ojalvo/ggHttMC-$1-SUBRct/dags/dag \
ggHttMC-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=False data=False

rm SUBRct.py

