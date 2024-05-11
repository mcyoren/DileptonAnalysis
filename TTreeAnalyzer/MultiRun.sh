#!/bin/sh

tcsh -c "source source.csh"

LIST=`ls -lhtr /phenix/plhf/mitran/taxi/Run14AuAu200CAMBPro109/19294/data/my-4*.root | awk '{printf("%s\n",$9)}'`

NUM=0
Additional_Rejection=0
Gen_Cut=0
isERT=1
#export MYINSTALL=/direct/phenix+u/vdoomra/install
#export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR

echo $DIR

for file in $LIST
do
  if (( $NUM == $1 ))
  then

    NAME=`echo $file | awk -F \/ '{printf("%s\n",$7)}'`
    NAME0=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    DUMMY=`echo $file | awk -F \/ '{printf("%s\n",$7)}' | awk -F\- '{printf("%s\n",$2)}' | awk -F\. '{printf("%s\n",$1)}'`

    DIRECTORY=$(( 10#$DUMMY ))
    output=first_$NAME0
    echo $file $output $DIRECTORY

    root -l -b <<EOF
    .L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/offline/AnalysisTrain/Run14AuAuLeptonComby/MyEvent.C+
    .L /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/NewHitAssociation.C+
    NewHitAssociation("$file","$output")
EOF

mv my-$output /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/First
  fi
  NUM=$(( $NUM + 1 ))
done

popd
rm -r $DIR
