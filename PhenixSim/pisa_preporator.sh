#!/bin/sh

mypath=$PWD
Green='\033[0;32m'
Red='\033[0;31m' 
Color_Off='\033[0m'

ReWriteVTX=0
if test "$#" -ne 2; then
    echo -e "   ${Red}pisa_preparator.sh run_number run ${Color_Off}"
    echo -e "   ${Red}run_number = number of the run (choose any run group), e.g. 406541 ${Color_Off}"
    echo -e "   ${Red}run = year of the run, e.g. 14 ${Color_Off}"
    echo -e "   ${Red}if you need vtx - add path to vtx source at third parameter: (run pisa to get those files) ${Color_Off}"
    if test "$#" -ne 3; then
        exit -1
    else 
        PathToVTX=$3
        ReWriteVTX=1
    fi
fi
RUNNO=$1
run=$2

printf "${Green}Your runnumber is   ${RUNNO} ${Color_Off}\n"
printf "${Green}Your run year  is   ${run} ${Color_Off}\n"

if test -d make_sim; then
printf "${Green}make_sim dir exists${Color_Off}\n"
else 
printf "${Green}creating make_sim${Color_Off}\n"
mkdir make_sim
fi

if test -d make_sim/pisa; then
printf "${Green}make_sim/pisa dir exists${Color_Off}\n"
else 
printf "${Green}creating pisa${Color_Off}\n"
mkdir make_sim/pisa
fi


cd make_sim/pisa

printf "${Green}getting all PISA source files${Color_Off}\n"

LuxorLinker.pl -1 $RUNNO #doing nothing: evrything is rewritting by 3 lines below....
/afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run${run}/pisaLinker.csh
/afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run${run}/pisaToDSTLinker.csh
./pisaToDSTLinker.csh ##sometimes is not working.... but line below alway works

cp /afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run${run}/pisaToDST_200GeV.C pisaToDST.C

#lines below are intended for vtx: those files are obtain by running pisa, if you need those - please run pisa and copy maniually 

if test "$ReWriteVTX" -ne 0; then
    printf "${Green}rewriting vtx files${Color_Off}\n"
    cp $PathToVTX/pisaToDST_VTX.C .
    cp $PathToVTX/glogon.kumac .
    cp $PathToVTX/pisa.input .
    cp $PathToVTX/pisa.kumac .  ## use your run number!!!!!!
    cp $PathToVTX/pisaToDST_IOManager.C .
    cp $PathToVTX/svx* .
fi 

printf "${Green}jumping back to initial folder and copying usefull scripts from git${Color_Off}\n"
cd $mypath

if test -d DileptonAnalysis; then
    printf "${Green}DileptonAnalysis exists -  redownloading${Color_Off}\n"
    rm -rf DileptonAnalysis
fi

printf "${Green}copying everything from git${Color_Off}\n"

git clone -n --depth=1 --filter=tree:0 \
  https://github.com/mcyoren/DileptonAnalysis
cd DileptonAnalysis
git sparse-checkout set --no-cone PhenixSim
git checkout

cd $mypath

printf "${Green}copying all necessary files from cloned folder${Color_Off}\n"
gitsim=DileptonAnalysis/PhenixSim
mv $gitsim/*.* .
mv $gitsim/embed .
mv $gitsim/embedreco .
mv $gitsim/real .
mv $gitsim/sim .

printf "${Green}Done! Go work!${Color_Off}\n"
#mv DileptonAnalysis/PhenixSim/* .

#cp ../glogon.kumac .
#cp ../pisa.input .
#cp ../pisaToDST.C .
