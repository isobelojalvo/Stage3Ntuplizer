#!/bin/sh                                                                                                                                                                                                                                                                                                                     
voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat testL1TRCTRatesEff.py > SUBRct.py
cat submitSecondaryFiles.py >> SUBRct.py

rm -r /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/
mkdir /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/
#submit dir = /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/submit
#make dag dir
mkdir -p /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/dags
mkdir -p /data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/submit --output-dag-file=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2c2/server?SFN=/hdfs/store/user/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/ ZeroBias-560627   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=QCD_Pt-15to3000-MINIAOD.txt \
--submit-dir=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/submit \
--output-dag-file=/data/ojalvo/QCD_Pt-15to3000MC-Jets-5-SUBRct/dags/dag \
QCD_Pt-15to3000MC-Jets-5  \
--vsize-limit=7000 \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/test/SUBRct.py farmout=True rates=False data=False jets=True

rm SUBRct.py

