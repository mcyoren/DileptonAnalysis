#!/bin/sh

itype=$2  #0 -- only single, 1 -- only pisa, 2 --only pisaToDST, 3 -- all
shift=1

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs/
#:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/vborisov/Pasha/install
#export INSTALL=/phenix/plhf/vborisov/Pasha/install
#export LD_LIBRARY_PATH=$INSTALL/lib:$LD_LIBRARY_PATH
#export TSEARCHPATH=/phenix/plhf/vborisov/Pasha/install
#export DCACHE_DOOR=phnxcore03.rcf.bnl.gov:1094
#export GSEARCHPATH ${GSEARCHPATH}:DCACHE

#echo $LD_LIBRARY_PATH

if test -d /phenix/plhf/mitran/Simul/Dileptons/output_single; then
echo "output exist"
else 
mkdir output
chmod g+rx /phenix/plhf/mitran/Simul/Dileptons/output
fi 

cd /phenix/plhf/mitran/Simul/Dileptons/output

if test -d /phenix/plhf/mitran/Simul/Dileptons/output_single; then
echo "output_single dir exist"
else 
echo "creating exist"
mkdir /phenix/plhf/mitran/Simul/Dileptons/output_single
fi

RUNNUMBERS_RG0=(409471 409471 409471)
RANDOM=$$
run_number=${RUNNUMBERS_RG0[ ($RANDOM) % ${#RUNNUMBERS_RG0[@]} ]}

echo "=====Generating a Random Run Number=========="
echo $run_number

echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="

INPUT=$(( $1 + $shift ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR
echo start the single simulation in directory $DIR

pwd
NEVT=10000

if [[ $itype == 0 || $itype == 3 ]];then
cp /phenix/plhf/mitran/Simul/Dileptons/make_sim/make_single.C .
cp /phenix/plhf/mitran/Simul/Dileptons/make_sim/heau200_bbcz.root .
root -l -b -q make_single.C

NLINE=`wc -l < oscar.particles.dat`

echo $NEVT events generated

if test -d /phenix/plhf/mitran/Simul/Dileptons/output_single; then
echo "output_single dir exist"
else 
echo "creating output_single"
mkdir /phenix/plhf/mitran/Simul/Dileptons/output_single
fi 

if test -d /phenix/plhf/mitran/Simul/Dileptons/output_single/single; then
echo "single dir exist"
else 
echo "creating single"
mkdir /phenix/plhf/mitran/Simul/Dileptons/output_single/single/
fi 

cp oscar.particles.dat /phenix/plhf/mitran/Simul/Dileptons/output_single/single/$DIR.oscar.particles.dat


else 
echo "no single"
fi


# add proper header to oscar file

#cp /phenix/plhf/malik/simulation/make_sim/oscar_header.txt .
#cat oscar.particles.dat >> oscar_header.txt
#mv oscar_header.txt oscar.particles.dat

#if [ -a /phenix/plhf/roli/embed/taxi/TTree_out_single_$DIR.root ]
#then
#  echo "simulation file already generated"
#  exit
#fi

mkdir pisa
pushd pisa

if [[ $itype == 1 || $itype == 3 ]];then

echo "==============================================="
echo "============= START PISA NOW =================="
echo "==============================================="

# set up Run14 PISA
mv ../oscar.particles.dat .
cp /phenix/plhf/mitran/Simul/Dileptons/make_sim/pisa/* .
#cp /phenix/plhf/mitran/Simul/Dileptons/pisa_svx/* .
#run_number=428737

sed -i 's/ptrig [0-9]*/ptrig '$NEVT'/' glogon.kumac

pisa<pisa.input

else 
echo "no pisa"
fi


if [[ $itype == 2 || $itype == 3 ]];then
echo "==============================================="
echo "================ PISA TO DST =================="
echo "==============================================="
SMEAR=1.0 #BBC resolution 0-20,20-40,40-60,60-93,2.17 for 0-100

#cp /phenix/plhf/mitran/Simul/Dileptons/make_sim/pisa/pisaToDST_VTX.C .
#cp /phenix/plhf/mitran/Simul/Dileptons/make_sim/pisa/pisaToDST_IOManager.C .
#sed -i '111s/.*/  vtx_sim->ZVertexSigma( '$SMEAR' );/' pisaToDST.C 
#root -b -q -l pisaToDST.C\($NEVT\)
root -b -q -l 'pisaToDST_VTX.C('$NEVT', "PISAEvent.root", "dst_out.root", "svxeval.root", '$run_number', 0)'

if test -d /phenix/plhf/mitran/Simul/Dileptons/output_single/dst; then
echo "dst dir exist"
else 
echo "creating dst dir"
mkdir /phenix/plhf/mitran/Simul/Dileptons/output_single/dst/
fi 

if test -f dst_out.root; then
cp dst_out.root /phenix/plhf/mitran/Simul/Dileptons/output_single/dst/dst_out_single_$DIR.root

popd
rm -rf pisa

popd
rm -r $DIR

echo "Done! Go Work!"

else 
echo "smt wrong!"
fi

else 
echo "no pisaToDST"
fi
