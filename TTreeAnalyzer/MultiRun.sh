#!/bin/sh

LIST=`ls -lhtr /phenix/plhf/mitran/taxi/Run14AuAu200CAMBPro109/19???/data/my-4*.root | awk '{printf("%s\n",$9)}'`

NUM=0
Additional_Rejection=0
Gen_Cut=0
isERT=1

echo $LD_LIBRARY_PATH
echo $TSEARCHPATH

tmpdir="/home/tmp/"${USER}"_jobreal_"$1

if test -d $tmpdir; then
echo $tmpdir exists
else 
mkdir -p $tmpdir
fi
echo "cd $tmpdir"
cd       $tmpdir


INPUT=$(( $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR

echo $DIR

mystart=`date +%s`

for file in $LIST
do
  if (( $NUM == $1 ))
  then

    NAME=`echo $file | awk -F \/ '{printf("%s\n",$7)}'`
    NAME0=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    DUMMY=`echo $file | awk -F \/ '{printf("%s\n",$7)}' | awk -F\- '{printf("%s\n",$2)}' | awk -F\. '{printf("%s\n",$1)}'`

    DIRECTORY=$(( 10#$DUMMY ))
    output=first_"$NAME"_"$NAME0"
    echo $file $output $DIRECTORY

    root -l -b <<EOF
    gSystem->Load("libRun14AuAuLeptonEvent");    
    gSystem->Load("/phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/NewHitAssociation_C.so");
    NewHitAssociation("$file","$output")
EOF

mv my-$output /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/All
  fi
  NUM=$(( $NUM + 1 ))
done

popd
rm -r $DIR
rm -r $tmpdir

myend=`date +%s`

mytot=$((myend-mystart))
mymin=60

mysec=$((mytot%mymin))
mymins=$((mytot/mymin))
myhour=$((mymins/mymin))

echo "Accumulated time:" $myhour "hours" $mymins "minutes" $mysec "seconds"

#gSystem->Load("/phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/offline/AnalysisTrain/Run14AuAuLeptonComby/MyEvent_C.so");
