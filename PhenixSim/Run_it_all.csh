#!/usr/local/bin/tcsh -f


#Do the real, sim, and embedding

#########################################################################
# This is based on the original logic by Takashi.                       #
#                                                                       #
# The idea is to set up the environment in a uniform way to facilitate  #
# accessing data among the various stages of the full embedding process.#
#                                                                       #
# You can find all sources on git                                       #
# https://github.com/mcyoren/DileptonAnalysis                           #
# If you have any questions please contact me:                          #
# Yuri Mitrankov (mitrankovy@gmail.com)                                 #
#########################################################################


setenv prompt 1
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
source $i
end
source $HOME/.login
unsetenv OFFLINE_MAIN
unsetenv ONLINE_MAIN
unsetenv ROOTSYS

if( $#argv != 6) then
echo "   Run_it_all.csh num"
echo "   num       = job number"
echo "   type of sim = single, pisa, embed, all sim (no embed), all"
echo "   type of input  = single, helios, pythia8"
echo "   type of particle  = photon, other"
echo "   shift  = 1-10000 - Nev"
echo "   Number of events  = 1-10000"
exit -1
endif

set jobno = $1
set selected_paticle = $4
set shift = $5
set NEVT = $6

echo "jobno is set to     $jobno              "
echo "type of sim         $2                  "
echo "type of input       $3                  "
echo "selected_paticle    $selected_paticle   "
echo "shift               $shift              "
echo "NEVT                $NEVT               "

setenv DATADIR $PWD

source /opt/phenix/bin/phenix_setup.csh new
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH $DATADIR/embedreco/install/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH $DATADIR/embed/svx_cent_ana/install/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH .:/phenix/hhj/hachiya/15.08/embed/embed:$TSEARCHPATH
setenv TSEARCHPATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/:$TSEARCHPATH

echo "LD_LIBRARY_PATH: "
echo $LD_LIBRARY_PATH
echo "TSEARCHPATH: "
echo $TSEARCHPATH

set vertex_txt_dir = $DATADIR/real/work/output/vertexes.txt
set outputsingle_dir = $DATADIR/output_single

set INPUT = `expr $shift + $jobno`
echo "input is " $INPUT

set DIR = `printf "%05d" $INPUT`

echo "INPUT               $INPUT              "
echo "DIR                 $DIR                "
echo "vertex_txt_dir      $vertex_txt_dir     "
echo "outputsingle_dir    $outputsingle_dir   "

if( -d $outputsingle_dir ) then
echo "outputsingle_dir exists"
else 
echo "mkdir $outputsingle_dir"
mkdir $outputsingle_dir
chmod g+rx $outputsingle_dir
endif

set RUNNUMBERS_RG0 = (409471 409471 409471)
set RANDOM=$$
set rnd = `expr $RANDOM % ${#RUNNUMBERS_RG0} + 1` ### for bash use like this and numerign from 0 ${#RUNNUMBERS_RG0[@]}
set run_number = $RUNNUMBERS_RG0[$rnd]

echo "=======Generating a Random Run Number=========="
echo "                " $run_number

echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="

set inputvtx = $DATADIR/real/work/output/vertexes.txt
set oscarname = $DIR.oscar.particles.dat

if( ( $2 == 0 || $2 == 3 || $2 == 4 ) && $3 == 0 ) then
echo "==============================================="
echo "============= SINGLE SIM STARTS ==============="
echo "==============================================="
set scriptdir   = $DATADIR/sim/gen/single
set macroname   = make_single.C
set outsingle   = $DATADIR/output_single/single
set tmpdir      = "/home/tmp/${USER}_job_$INPUT"
set ptmin = 0.4
set ptmax = 10.0
set n     = -1. #n: <0 hagdorn (mb HeAu), =0 flat, >0 power law
set id    = $selected_paticle #0,1,2,3,4-pi0,pi+,pi-,e+,e-################helios jpsi and phi####pythia ccbar bbar

echo "jobno          $jobno            "
echo "run_number     $run_number       "
echo "ptmin          $ptmin            "
echo "ptmax          $ptmax            "
echo "n              $n                "
echo "id             $id               "
echo "inputvtx       $inputvtx         "
echo "oscarname      $oscarname        "
echo "scriptdir      $scriptdir        "
echo "macroname      $macroname        "
echo "outsingle      $outsingle        "
echo "tmpdir         $tmpdir           "

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir

cp $scriptdir/$macroname .
echo "start running root -l -b -q" 'make_single.C("'$vertex_txt_dir'","'$oscarname'",'$NEVT','$ptmin','$ptmax','$n','$id')'
root -l -b -q 'make_single.C("'$vertex_txt_dir'","'$oscarname'",'10000','$ptmin','$ptmax','$n','$id')'
echo "finished rinnung"
echo "copying"
cp $oscarname  $outsingle/
echo "finished copying"
#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "removed $tmpdir"

endif

if( ( $2 == 0 || $2 == 3 || $2 == 4 ) && $3 == 1 ) then
echo "==============================================="
echo "============= HELIOS TO OSCAR ================="
echo "==============================================="
set inputhelios = $DATADIR/output_single/helios/helios_phi_ee_02_5_10M.root
if( $selected_paticle == 2) then
 set inputhelios = $DATADIR/output_single/helios/helios_jpsi_ee_02_5_11M.root
endif
if( $selected_paticle == 0) then
 set inputhelios = $DATADIR/output_single/helios/helios_pi0_gg_3_10_50M.root
endif
if( $selected_paticle == 3) then
 set inputhelios = $DATADIR/output_single/helios/helios_pi0_gee_05_2_10M.root
endif
set scriptdir   = $DATADIR/sim/gen/HELIOS/work
set scriptname  = Convert_HELIOS.csh
set macroname   = WriteROOT2Oscar.C
set outsingle   = $DATADIR/output_single/helios
set tmpdir      = "/home/tmp/${USER}_job_$INPUT"

echo "jobno          $jobno            "
echo "run_number     $run_number       "
echo "inputhelios    $inputhelios      "
echo "inputvtx       $inputvtx         "
echo "oscarname      $oscarname        "
echo "scriptdir      $scriptdir        "
echo "scriptname     $scriptname       "
echo "macroname      $macroname        "
echo "outsingle      $outsingle        "
echo "tmpdir         $tmpdir           "

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir

cp $scriptdir/$scriptname .
cp $scriptdir/$macroname .
echo running "./$scriptname $inputvtx $jobno 10000 $inputhelios $oscarname"
./$scriptname $inputvtx $jobno 10000 $inputhelios $oscarname
echo "finished rinnung"
echo "copying"
cp $oscarname  $outsingle/
echo "finished copying"
#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "removed $tmpdir"

endif

if( ( $2 == 0 || $2 == 3 || $2 == 4 ) && $3 == 2 ) then
echo "==============================================="
echo "============= PYTHIA TO OSCAR ================="
echo "==============================================="
set inputpythia = $DATADIR/output_single/pythia8/ccbartree$DIR.root
if( $selected_paticle == 2) then
 set inputpythia = $DATADIR/output_single/pythia8/bbbartree$DIR.root
endif
set scriptdir   = $DATADIR/sim/gen/pythia8
set scriptname  = Convert_pythia8.csh
set macroname   = WriteROOT2OscarPythia.C
set outsingle   = $DATADIR/output_single/pythia8
set tmpdir      = "/home/tmp/${USER}_job_$INPUT"

echo "jobno          $jobno            "
echo "run_number     $run_number       "
echo "inputpythia    $inputpythia      "
echo "inputvtx       $inputvtx         "
echo "oscarname      $oscarname        "
echo "scriptdir      $scriptdir        "
echo "scriptname     $scriptname       "
echo "macroname      $macroname        "
echo "outsingle      $outsingle        "
echo "tmpdir         $tmpdir           "

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir

cp $scriptdir/$scriptname .
cp $scriptdir/$macroname .
echo running "./$scriptname $inputvtx $inputpythia $oscarname"
./$scriptname $inputvtx $inputpythia $oscarname
echo "finished rinnung"
echo "copying"
cp $oscarname  $outsingle/
echo "finished copying"
#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "removed $tmpdir"

endif


if( $2 == 1 || $2 == 3 || $2 == 4 ) then

echo "==============================================="
echo "============= START PISA NOW =================="
echo "==============================================="

if ( $3 == 1 ) then 
set oscarname = $DATADIR/output_single/helios/$DIR.oscar.particles.dat
else if ( $3 == 2 ) then 
set oscarname = $DATADIR/output_single/pythia8/$DIR.oscar.particles.dat
else 
set oscarname = $DATADIR/output_single/single/$DIR.oscar.particles.dat
endif

set pisadir  = $DATADIR/make_sim/pisa

set outpisadir  = $DATADIR/output_single/simdst
set outdstname  = dst_out_single_$DIR.root
set tmpdir      = "/home/tmp/${USER}_pisa_$INPUT"

echo "jobno          $jobno            "
echo "run_number     $run_number       "
echo "Nev            $NEVT             "
echo "oscarname      $oscarname        "
echo "pisadir        $pisadir          "
echo "outpisadir     $outpisadir       "
echo "tmpdir         $tmpdir           "

#checkng outdir
if( ! -d $outpisadir ) then
echo "creating $outpisadir" 
mkdir -p $outpisadir
endif
#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir


cp $oscarname oscar.particles.dat
cp $pisadir/* .

sed -i 's/ptrig [0-9]*/ptrig '$NEVT'/' glogon.kumac

echo running "pisa<pisa.input"
pisa<pisa.input
echo "finished rinnung"

echo "==============================================="
echo "================ PISA TO DST =================="
echo "==============================================="

echo "running root -b -q -l pisaToDST_VTX.C"
root -b -q -l 'pisaToDST_VTX.C('$NEVT', "PISAEvent.root", "dst_out.root", "svxeval.root", '$run_number', 0)'
echo "finished rinnung"

if ( -f dst_out.root ) then
echo "moving"
mv dst_out.root $outpisadir/$outdstname
echo "finished moving"
else 
echo "something is fkng wrong"
endif

#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "removed $tmpdir"

endif

##########################################
################EMBEDDING#################
#################STARTS###################
##########################################
if( $2 == 2 || $2 == 4 ) then
echo "==============================================="
echo "================= EMBEDDING ==================="
echo "==============================================="

setenv sourcedir /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/embed/svx_cent_ana/build
if( $1 == 0) then
cd $sourcedir
make -j4
make install
endif

echo "Setting internal variables..."

set runnum    = $run_number
set jobno     = $1
set evtnum    = $NEVT

set INPUT = `expr $shift + $jobno`
echo $INPUT

set outntana  = tree"$DIR".root

if ( $3 == 1 ) then 
set inputoscar = $DATADIR/output_single/helios/$DIR.oscar.particles.dat
else if ( $3 == 2 ) then 
set inputoscar = $DATADIR/output_single/pythia8/$DIR.oscar.particles.dat
else 
set inputoscar = $DATADIR/output_single/single/$DIR.oscar.particles.dat
endif
set inputsim = $DATADIR/output_single/simdst/dst_out_single_$DIR.root
set inputvtx = $DATADIR/real/work/output/vertexes.txt
set inputreal = $DATADIR/real/work/outputnew/CNTmerge_MB-0000$run_number-0000.root
set outdst    = kek0.root
set embedpartID = $selected_paticle

set scriptdir = $DATADIR/embed/work

set outmytreedir   = $DATADIR/output_single/embed

#checkng outdir
if( ! -d $outmytreedir ) then
echo "creating $outmytreedir" 
mkdir -p $outmytreedir
endif

set tmpdir      = "/home/tmp/${USER}_embed_$INPUT"

set inreal = `basename $inputreal`
set insim  = `basename $inputsim`

set runnum = `echo $inreal | awk -F'-' '{printf "%d", $2}'`



echo "runnum       $runnum         "
echo "jobno        $jobno          "
echo "evtnum       $evtnum         "
echo "inputreal    $inputreal      "
echo "inputsim     $inputsim       "
echo "outdst       $outdst         "
echo "embedpartID  $embedpartID    "
echo "outntana     $outntana       "
echo "inputvtx     $inputvtx       "
echo "inputoscar   $inputoscar     "
echo "scriptdir    $scriptdir      "
echo "outmytreedir $outmytreedir   "
echo "tmpdir       $tmpdir         "

# exit


#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir


##pisaToDSTLinker with new library
/opt/phenix/core/bin/LuxorLinker.pl -1 $runnum

cp ${scriptdir}/Fun4All_embedeval.C .
cp ${scriptdir}/Fun4All_embedeval_svx.C .
cp ${scriptdir}/embed_IOManager.C    .
cp ${scriptdir}/svxPISA.par         .
cp ${scriptdir}/svx_threshold.dat   .

#copy input file
cp $inputreal .
cp $inputsim .

echo "pwd"
pwd
ls -ltr

echo "yolo"
echo ".x Fun4All_embedeval_svx.C($evtnum, "'"'$insim'"'", "'"'$inreal'"'", "'"'$outdst'"'", "$embedpartID", "'"'$outntana'"'", "'"'$inputvtx'"'", "'"'$inputoscar'"'", $runnum);" >  cmd.input
echo ".q" >> cmd.input

##run root
root -b < cmd.input #>& $HOME/root.log
#ls -ltr
#move
#mv -f $outdst    $outdstdir
#mv -f kek2.root $outmytreedir/
#mv -f $outntana  $outmytreedir/
mv -f my-$outntana $outmytreedir/

#remove tmp dir
cd $HOME
rm -fr $tmpdir
echo "removed $tmpdir"

endif
##########################################
###################END####################
################EMBEDDING#################
##########################################