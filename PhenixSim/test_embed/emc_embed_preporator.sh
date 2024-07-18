#!/bin/sh
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs/

itype=$1  #0 -- all; 1 -- create realdata list, 2 -- split data, 3 -delete leftovers

chmod g+rx /phenix/plhf/mitran/Simul/ETA/output
cd /phenix/plhf/mitran/Simul/ETA/output

echo "==============================================="
echo "============= START PREPARATION ==============="
echo "==============================================="


if [[ $itype == 1 || $itype == 0 ]];then

echo "==============================================="
echo "=========== Creating realdata list ============"
echo "==============================================="

ls -Sr /pnfs/rcf.bnl.gov/phenix/phnxreco/run14/run14HeAu_200GeV_CA_pro102/run_0000416000_0000417000/CNT_MB/CNT_MB_run14HeAu_200GeV_CA_pro102* > realdata_run_embed_list.txt

else 
echo "list already exist"
fi

CNTin=`cat realdata_run_embed_list.txt`

if test -d /phenix/plhf/mitran/Simul/ETA/output/realdatadst; then
echo "realdatadst dir exist"
else 
echo "creating exist"
mkdir realdatadst
fi

cd realdatadst
NUM=0
CntFileNames=()
for icnt in $CNTin
do
  echo $icnt
  if [[ $NUM -gt 9 ]]
  then
    echo "Files are downaloaded"
    break
  else
    iCNT=`ls $icnt | awk -F\/ '{printf("%s\n",$10)}'`
    CntFileNames+=($iCNT)
    echo $iCNT
    if test -f $iCNT; then
      echo "Great, $iCNT is in realdatadst"
    else 
      echo xrdcp root://phnxcore03.rcf.bnl.gov:1094$icnt .
      #xrdcp root://phnxcore03.rcf.bnl.gov:1094$icnt .
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

ls -Sshr . 
cd ..

if test -d /phenix/plhf/mitran/Simul/ETA/output_pi0/realdstfiles; then
echo "embed dir exist"
else 
echo "creating realdstfiles dir"
mkdir /phenix/plhf/mitran/Simul/ETA/output_pi0/realdstfiles
fi 

pushd /phenix/plhf/mitran/Simul/ETA/output_pi0/realdstfiles

echo start splitting proccess in /phenix/plhf/mitran/Simul/ETA/output_pi0/realdstfiles

if [[ $itype == 2 || $itype == 0 ]];then

echo "==============================================="
echo "============ START Splitting data ============="
echo "==============================================="

cp /phenix/plhf/mitran/Simul/ETA/embedding_setup_HeAu200/split_realdata.C .

ii=0
for value in "${CntFileNames[@]}"
do
     echo split $value to different centralities
     root -b -q 'split_realdata.C("'/phenix/plhf/mitran/Simul/ETA/output/realdatadst/$value'",'$ii', 100000,5,50000)'
      ii=$(( $ii + 1 ))
done

else 
echo "data arent splitted"
fi

#ls -Sr /phenix/plhf/mitran/Simul/ETA/output/realdatadst/* > realdata_run_embed_list.txt

if [[ $itype == 3 || $itype == 0 ]];then
echo "deleting dst data"
rm /phenix/plhf/mitran/Simul/ETA/output/realdatadst/*
else 
echo "dst data left"
fi