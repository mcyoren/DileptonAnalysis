#!/bin/tcsh

source /opt/phenix/bin/phenix_setup.csh
#setenv TSEARCHPATH /direct/phenix+u/workarea/ahodges/run14HadEff/install
setenv LD_LIBRARY_PATH /direct/phenix+u/ahodges/cvs/install/lib:${LD_LIBRARY_PATH}

#Inputs into the shell script
# $1 - process number
# $2 - number of event being processed in each job
echo "input parameters--- "
echo "process number:   " $1
echo "number of events: " $2
echo "particle: " $3


set runPISA=0       # 1: run PISA and everything, 0: run reconstruction and ana code only
set doPtoDST=0        # 1: re-analyze PISAEvent.root if there were changes to PISAtoDST.C
#----------------------------------------------------------------#
# go to my simulation directory, then copy pisa scripts over
#----------------------------------------------------------------# 
#simulation stuff store here:
set simDir=/gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff/simOutput14
cd $simDir
pwd

if( ! -d job ) then
    mkdir job
endif

if( ! -d job/job_$3_$1 ) then
    mkdir job/job_$3_$1
endif

#----------------------------------------------------------------#
# run PISA
#----------------------------------------------------------------#
cd job/job_$3_$1

if ( $runPISA ) then
   #cp -r /gpfs/mnt/gpfs02/phenix/hhj/hhj1/cwong14/run15hadronEff/pisa .   #!!!!! skip this line if not running PISA !!!!
    cp -r /gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff/pisaStuff/run14 .
    cd run14/ 
endif

  
#cp /gpfs/mnt/gpfs02/phenix/hhj/hhj1/cwong14/run15hadronEff/pisa/pisaToDST_200GeV_run15Deadmap_ISU.C .
#cp /gpfs/mnt/gpfs02/phenix/hhj/hhj1/cwong14/run15hadronEff/pisa/DchEfficiency0428267_combined.Real .
#----------------------------------------------------------------#
# modify glogon.kumac
# need to find number of events that this job should skip s.t. each job is unique
# Process number * number of events in each job
#----------------------------------------------------------------#
if ( $runPISA ) then
    @ nskip = $1 * $2
    echo "===================skipping" $nskip "entries====================="
    
    #ln -sf /direct/phenix+u/cwong14/run10_run11_Correlation/run15pp/hadronEff/oscarInput/$3.root oscar.root
    ln -sf /gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/oscarFiles/run14/$3.root oscar.root
    echo "oscar oscar.root 1">glogon.kumac
    echo "nskip "$nskip"">>glogon.kumac
    echo "ptrig "$2>>glogon.kumac
    echo "exit">>glogon.kumac
    echo "================glogon.kumac written for oscarroot==============="
    
    echo "=========================== run PISA ============================"
    pisa < pisa.input 
endif

if ( $doPtoDST ) then
#convert pisa output to DST format so Fun4All can read it
echo "=========== convert PISAEvent.root to simDST files =============="
#root -l -b -q pisaToDST_200GeV_run7Deadmap.C
#root -l -b -q pisaToDST_200GeV_DCHdeadmap.C       #with run12 DCH deadmap
#root -l -b -q pisaToDST_200GeV.C                  #with run 9 DCH deadmap
#root -l -b -q pisaToDST_200GeV_run12Deadmap_ISU.C
#root -l -b -q pisaToDST_200GeV_run15Deadmap_ISU.C
echo $PWD
#root -l -b -q 'pisaToDST_200GeV.C('$2','\"'PISAEvent.root'\"','\"'simDST.root'\"',405863)'
root -l -b -q 'pisaToDST_200GeV.C('$2','\"'PISAEvent.root'\"','\"'simDST.root'\"',414507)'
mv simDST.root $simDir/DST_rootFile/simDST_$3_$1.root
endif
#----------------------------------------------------------------#
# make QA histo
#----------------------------------------------------------------#
echo "========================= make histo ========================="
cd $simDir
#root -l -b -q /gpfs/mnt/gpfs02/phenix/hhj/hhj1/cwong14/run15hadronEff/ana/src/RunThisMacro.C'('\"'DST_rootFile/simDST_'$3'_'$1'.root'\"','\"'histoQA_'$3'_'$1'.root'\"','\"'../pisa/Run15ppDataPhi.root'\"')'
root -l -b -q /gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff/analysis/RunThisMacro.C'('\"'DST_rootFile/simDST_'$3'_'$1'.root'\"','\"'histoQA_'$3'_'$1'.root'\"',14,1,0)'
mv histoQA_$3_$1.root histoFile/histoQA_$3_$1.root

#----------------------------------------------------------------#
# remove pisa scripts, comment out for debugging one run, otherwise
# uncomment to save disk space
#----------------------------------------------------------------#
    

rm -r job/job_$3_$1

