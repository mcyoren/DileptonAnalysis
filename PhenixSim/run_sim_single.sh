#!/bin/sh

itype=$2  #0 -- only single, 1 -- only pisa, 2 --only pisaToDST, 3 -- all
shift=0

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

RUNNUMBERS_RG0=(405863 406268 406773 407176 407608 407951 405864 406539 406774 407195 407610 407953 405865 406540 406831 407196 407611 407959 405868 406543 406832 407197 407614 407960 405869 406544 406833 407198 407618 407963 405961 406546 406848 407199 407620 407964 405962 406549 406849 407200 407621 407965 405964 406555 406850 407269 407660 407966 405966 406579 406851 407270 407661 408070 405971 406581 406853 407271 407662 408071 405972 406584 406857 407272 407664 408073 405973 406661 406858 407276 407666 408074 405975 406662 406859 407362 407669 408075 405977 406663 406861 407367 407670 408076 405981 406666 406866 407368 407671 408077 405982 406671 406867 407369 407672 408078 405984 406674 406869 407370 407673 408175 405996 406675 406870 407371 407676 408176 406087 406676 406871 407372 407711 408177 406088 406677 406872 407375 407712 408181 406089 406697 406881 407376 407786 408182 406090 406698 406882 407377 407790 408183 406093 406700 406887 407378 407792 408184 406095 406745 406889 407379 407796 408217 406099 406746 406891 407380 407797 408218 406180 406747 406893 407381 407798 408219 406181 406753 406898 407445 407799 408220 406182 406754 406902 407447 407800 408224 406183 406759 406905 407448 407802 408225 406190 406760 406920 407454 407839 408226 406257 406761 406921 407455 407842 408227 406258 406762 407143 407456 407944 408229 406259 406764 407144 407457 407945 406261 406766 407145 407523 407946 406263 406767 407146 407524 407947 406265 406769 407147 407526 407948 406266 406772 407175 407558 407950)
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
NEVT=100

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
cp /phenix/plhf/mitran/Simul/Dileptons/pisa_svx/* .
run_number=428737

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
#rm -rf pisa

popd
#rm -r $DIR

echo "Done! Go Work!"

else 
echo "smt wrong!"
fi

else 
echo "no pisaToDST"
fi
