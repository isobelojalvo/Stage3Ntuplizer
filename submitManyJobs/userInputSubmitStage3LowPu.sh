#!/bin/sh                                                                                                                              

#voms-proxy-init --voms cms --valid 100:00                                                                                              

cat $1 > SUBStage3.py
cat submit.py >> SUBStage3.py

rm -r /data/ojalvo/ZeroBias-275782-$2-SUBStage3/
mkdir /data/ojalvo/ZeroBias-275782-$2-SUBStage3/
##make dag dir
mkdir -p /data/ojalvo/ZeroBias-275782-$2-SUBStage3/dags/daginputs

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=Run2016C-275782-ZeroBias-RAW.txt \
--submit-dir=/data/ojalvo/ZeroBias-275782-$2-SUBStage3/submit \
--output-dag-file=/data/ojalvo/ZeroBias-275782-$2-SUBStage3/dags/dag \
ZeroBias-275782-$2  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/submitManyJobs/SUBStage3.py farmout=True rates=True data=True $3 $4

rm SUBStage3.py

