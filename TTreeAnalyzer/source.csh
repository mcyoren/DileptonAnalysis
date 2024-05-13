#!/bin/csh

source /opt/phenix/bin/phenix_setup.csh -n new

setenv LD_LIBRARY_PATH .:/opt/phenix/bin:$LD_LIBRARY_PATH

#setenv LD_LIBRARY_PATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/lib:$LD_LIBRARY_PATH
#setenv TSEARCHPATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/:$TSEARCHPATH

#echo $LD_LIBRARY_PATH
#echo $TSEARCHPATH

#root -l -b <<EOF
#.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/offline/AnalysisTrain/Run14AuAuLeptonComby/MyEvent.C+
#.L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/NewHitAssociation.C+
#EOF

