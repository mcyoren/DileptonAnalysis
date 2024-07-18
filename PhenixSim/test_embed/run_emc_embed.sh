#!/bin/sh

shift=400
INPUT=$(( $1 + $shift ))
iSIM=`printf "%05d" $INPUT`
iCNT=$(( $1 % 10 ))
itype=$2  #0 -- run embed, 1 -- run embed without centrality


chmod g+rx /phenix/plhf/mitran/Simul/ETA/output
cd /phenix/plhf/mitran/Simul/ETA/output

echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="


# create new dir
DIR=`printf "%06d" $INPUT`
echo $DIR

if test -d /phenix/plhf/mitran/Simul/ETA/output/$DIR; then
echo "embed dir exist"
else 
echo "creating dir"
mkdir $DIR
fi 

if test -d /phenix/plhf/mitran/Simul/ETA/output_pi0/embed; then
echo "embed dir exist"
else 
echo "creating dir"
mkdir /phenix/plhf/mitran/Simul/ETA/output_pi0/embed
fi 

pushd $DIR

echo start the emdedding proccess in $DIR

export SIMFile=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/ETA/output_pi0/dst/dst_out_pi0_$iSIM.root
export CntFile=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/ETA/output_pi0/realdstfiles/realdst$iCNT
echo $CntFile
cp $CntFile* . 
cp /phenix/plhf/mitran/Simul/ETA/embedding_setup_HeAu200/run_embed.C .
ln -s /phenix/plhf/mitran/Simul/ETA/embedding_setup_HeAu200/lookup_3D_one_phi.root .

if [[ $itype == 1 || $itype == 0 ]];then

echo "==============================================="
echo "============= START EMBEDDING NOW ============="
echo "==============================================="

for icent in {0..3}
do
	
  realdst='realdst'$iCNT'_'$icent'.root'
  #realdst='/phenix/plhf/mitran/Simul/ETA/output/realdatadst/CNT_MB_run14HeAu_200GeV_CA_pro102-0000416021-9008.root'
    
  ls -lrt
  
  echo generating embeded file embed_pi0-0-20GeV-$icent'_'$iSIM.root with $realdst and $SIMFile
  #echo 'run_embed.C("'$realdst'","'$SIMFile'","embed_pi0-0-20GeV-'$icent'_'$INPUT'.root", 0, 10000)'
  root -b -q 'run_embed.C("'$realdst'","'$SIMFile'","embed_pi0-0-20GeV-'$icent'_'$iSIM'.root", 20000, 50000)'

  if test -f embed_pi0-0-20GeV-$icent'_'$iSIM.root; then
    echo "Great, go work!"
    mv embed_pi0-0-20GeV-$icent'_'$iSIM.root /phenix/plhf/mitran/Simul/ETA/output_pi0/embed
    if [[ $icent -gt 2 ]]
    then
      popd
      rm -r $DIR
    else 
      echo FAIL
    fi
  else 
    echo "finised in" 
  fi

done

else 
echo "embeddeding delayed!"
fi