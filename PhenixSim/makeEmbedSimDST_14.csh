#!/bin/tcsh

#source /opt/phenix/bin/phenix_setup.csh
#setenv TSEARCHPATH /direct/phenix+u/workarea/ahodges/run14HadEff/install
#setenv LD_LIBRARY_PATH /direct/phenix+u/ahodges/cvs/install/lib:${LD_LIBRARY_PATH}

#Inputs into the shell script
# $1 - process number
# $2 - number of event being processed in each job
echo "input parameters--- "
echo "process number:   " $1
echo "number of events: " $2
echo "particle: " $3
echo "run number: " $4

@ nskip = $1 * $2
echo "===================skipping" $nskip "entries====================="

cp -r /gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff/pisaStuff/run$4/* .


ln -sf /gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/oscarFiles/$3.root oscar.root
echo "oscar oscar.root 1">glogon.kumac
echo "nskip "$nskip"">>glogon.kumac
echo "ptrig "$2>>glogon.kumac
echo "exit">>glogon.kumac
echo "================glogon.kumac written for oscarroot==============="
    
echo "=========================== run PISA ============================"
pisa < pisa.input 

root -l -b -q 'pisaToDST_200GeV.C('$2','\"'PISAEvent.root'\"','\"'simDST.root'\"',414507)'

mv simDST.root simDST_$3_$1_$2_Events.root
