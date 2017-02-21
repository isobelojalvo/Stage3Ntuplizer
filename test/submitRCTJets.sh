#!/bin/sh                                                                                                                              
voms-proxy-init --voms cms --valid 100:00                                                                                              

cat testL1TRCTRatesEff.py > SUBRct.py
cat submit.py >> SUBRct.py

rm -r /data/ojalvo/ZeroBias-260627-Jets-SUBRct/
mkdir /data/ojalvo/ZeroBias-260627-Jets-SUBRct/
#submit dir = /data/ojalvo/ZeroBias-260627-Jets-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-260627-Jets-SUBRct/dags
mkdir -p /data/ojalvo/ZeroBias-260627-Jets-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-Jets-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-260627-Jets-SUBRct/submit --output-dag-file=/data/ojalvo/ZeroBias-260627-Jets-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-260627-Jets-SUBRct/ ZeroBias-260627-Jets-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-260627-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-260627-Jets-SUBRct/submit \
--output-dag-file=/data/ojalvo/ZeroBias-260627-Jets-SUBRct/dags/dag \
ZeroBias-260627-Jets  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=True data=True jets=True

rm -r /data/ojalvo/ZeroBias-259721-Jets-SUBRct/
mkdir /data/ojalvo/ZeroBias-259721-Jets-SUBRct/
#submit dir = /data/ojalvo/ZeroBias-259721-Jets-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/ZeroBias-259721-Jets-SUBRct/dags
mkdir -p /data/ojalvo/ZeroBias-259721-Jets-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-259721-Jets-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/ZeroBias-259721-Jets-SUBRct/submit --output-dag-file=/data/ojalvo/ZeroBias-259721-Jets-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/ZeroBias-259721-Jets-SUBRct/ ZeroBias-259721-Jets-v2c2   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2015D-259721-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-259721-Jets-SUBRct/submit \
--output-dag-file=/data/ojalvo/ZeroBias-259721-Jets-SUBRct/dags/dag \
ZeroBias-259721-Jets  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=True data=True jets=True

rm SUBRct.py

