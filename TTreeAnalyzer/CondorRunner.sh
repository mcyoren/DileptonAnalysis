#!/bin/sh

if [[ $1 == 0 ]];then
    rm /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/All/*
    rm /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/logs/*
fi

tcsh -c "source source.csh"

export LD_LIBRARY_PATH=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/lib:$LD_LIBRARY_PATH
export TSEARCHPATH=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install:$TSEARCHPATH

echo $LD_LIBRARY_PATH
echo $TSEARCHPATH

cd /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/build/
make -j4
make install

cd /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/ 

root -l -b <<EOF
gSystem->Load("libRun14AuAuLeptonEvent");
.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/NewHitAssociation.C+
EOF

if [[ $1 -gt 0 ]];then
    ./MultiRun.sh $1
fi

#.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/offline/AnalysisTrain/Run14AuAuLeptonComby/MyEvent.C+

if [[ $1 == 0 ]];then
    condor_submit multirun.job
fi