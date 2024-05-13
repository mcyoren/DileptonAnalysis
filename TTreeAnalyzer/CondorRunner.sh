#!/bin/sh

rm /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/All/*
rm /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/logs/*

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

#./MultiRun.sh 100
#.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/offline/AnalysisTrain/Run14AuAuLeptonComby/MyEvent.C+

condor_submit multirun.job
