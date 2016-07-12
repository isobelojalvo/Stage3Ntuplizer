#!/bin/sh                                                                                                                                                                                                                                                                                                                     
#voms-proxy-init --voms cms --valid 100:00
echo $1 $2 $3 $4
cat $1 > SUBStage3.py
cat submitSecondaryFiles.py >> SUBStage3.py

rm -r /data/ojalvo/ggHttMC-$2-SUBStage3/
mkdir /data/ojalvo/ggHttMC-$2-SUBStage3/
#make dag dir
mkdir -p /data/ojalvo/ggHttMC-$2-SUBStage3/dags/daginputs

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=ggHtt-76X-TSG-MINI.txt \
--submit-dir=/data/ojalvo/ggHttMC-$2-SUBStage3/submit \
--output-dag-file=/data/ojalvo/ggHttMC-$2-SUBStage3/dags/dag \
ggHttMC-$2  \
$CMSSW_BASE \
$CMSSW_BASE/src/L1Trigger/Stage3Ntuplizer/submitManyJobs/SUBStage3.py farmout=True rates=False data=False $3 $4 

rm SUBStage3.py

