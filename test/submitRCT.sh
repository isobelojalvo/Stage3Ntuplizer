#!/bin/sh                                                                                                                              
voms-proxy-init --voms cms --valid 100:00                                                                                              

cat testL1TRCTRatesEff.py > SUBRct.py
cat submit.py >> SUBRct.py

#Run2016C-276870-ZeroBias-RAW.txt

rm -r /data/ojalvo/ZeroBias-276870-SUBRct/
mkdir /data/ojalvo/ZeroBias-276870-SUBRct/
#submit dir = /data/ojalvo/ZeroBias-276870-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-276870-SUBRct/dags
mkdir -p /data/ojalvo/ZeroBias-276870-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-276870-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-276870-SUBRct/submit --output-dag-file=/data/ojalvo/ZeroBias-276870-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-276870-SUBRct/ ZeroBias-276870-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2016C-276870-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-276870-SUBRct/submit \
--output-dag-file=/data/ojalvo/ZeroBias-276870-SUBRct/dags/dag \
ZeroBias-276870  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=True data=True

rm SUBRct.py
exit;

rm -r /data/ojalvo/ZeroBias-259721-SUBRct/
mkdir /data/ojalvo/ZeroBias-259721-SUBRct/
#submit dir = /data/ojalvo/ZeroBias-259721-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-259721-SUBRct/dags
mkdir -p /data/ojalvo/ZeroBias-259721-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-259721-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-259721-SUBRct/submit --output-dag-file=/data/ojalvo/ZeroBias-259721-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-259721-SUBRct/ ZeroBias-259721-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-259721-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-259721-SUBRct/submit \
--output-dag-file=/data/ojalvo/ZeroBias-259721-SUBRct/dags/dag \
ZeroBias-259721  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=True data=True



rm -r /data/ojalvo/ZeroBias-260627-SUBRct/
mkdir /data/ojalvo/ZeroBias-260627-SUBRct/
#submit dir = /data/ojalvo/ZeroBias-260627-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-260627-SUBRct/dags
mkdir -p /data/ojalvo/ZeroBias-260627-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-260627-SUBRct/submit --output-dag-file=/data/ojalvo/ZeroBias-260627-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-SUBRct/ ZeroBias-260627-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-260627-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-260627-SUBRct/submit \
--output-dag-file=/data/ojalvo/ZeroBias-260627-SUBRct/dags/dag \
ZeroBias-260627  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=True data=True

rm SUBRct.py

