#!/bin/sh

echo "Start the routine"

if [[ $1 == 0 ]];then
    echo "deleting prev root files and logs"
    rm /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/All/*
    rm /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/logs/*
fi

echo "sourcing"

tcsh -c "source source.csh"

echo "finish sourcing and start setting lib and search path"

export LD_LIBRARY_PATH=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/lib:$LD_LIBRARY_PATH
export TSEARCHPATH=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install:$TSEARCHPATH

echo $LD_LIBRARY_PATH
echo $TSEARCHPATH

echo "start making myevent lib"

cd /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/build/
make -j4
make install

echo "Done making myevent lib"

cd /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/ 

echo "Begin compiling my script lib: .L NewHitAssociation.C+"

root -l -b <<EOF
gSystem->Load("libRun14AuAuLeptonEvent");
.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/NewHitAssociation.C+
EOF

echo "Done compiling"

if [[ $1 -gt 0 ]];then
    echo "runing one job:" $1      
    ./MultiRun.sh $1
fi

#.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/offline/AnalysisTrain/Run14AuAuLeptonComby/MyEvent.C+

if [[ $1 == 0 ]];then
    condor_rm -all
    echo "subbmiting 1000 jobs to condor"
    condor_submit multirun.job
fi

echo "DONE!!!"