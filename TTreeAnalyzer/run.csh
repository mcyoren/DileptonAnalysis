#!/bin/csh

source /opt/phenix/bin/phenix_setup.csh -n new

setenv LD_LIBRARY_PATH .:/opt/phenix/bin:$LD_LIBRARY_PATH

setenv LD_LIBRARY_PATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/:$TSEARCHPATH

echo $LD_LIBRARY_PATH
echo $TSEARCHPATH

set mystart=`date +%s`

root -l -b << EOF
    .L ../AnaTrain/Run14AuAuLeptonComby/MyEvent.C+
    .x NewHitAssociation.C+
EOF

set myend=`date +%s`

@ mytot=$myend - $mystart
set mymin=60

@ mysec=$mytot % $mymin
@ mymins=$mytot / $mymin
@ myhour=$mymins / $mymin

echo "Accumulated time:" $myhour "hours" $mymins "minutes" $mysec "seconds"