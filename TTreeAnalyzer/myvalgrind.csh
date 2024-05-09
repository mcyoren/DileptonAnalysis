#!/bin/csh

source /opt/phenix/bin/phenix_setup.csh -n new

setenv LD_LIBRARY_PATH .:/opt/phenix/bin:$LD_LIBRARY_PATH

setenv LD_LIBRARY_PATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/:$TSEARCHPATH

echo $LD_LIBRARY_PATH
echo $TSEARCHPATH

valgrind -v --num-callers=40 --leak-check=full --error-limit=no \
    --log-file=valgrind.log --suppressions=$ROOTSYS/root.supp \
    --leak-resolution=high root.exe

