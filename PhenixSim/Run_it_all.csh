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

setenv HOME /phenix/u/$LOGNAME
setenv prompt 1
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
source $i
end
source $HOME/.login
unsetenv OFFLINE_MAIN
unsetenv ONLINE_MAIN
unsetenv ROOTSYS

set Green='\033[0;32m'
set Red='\033[0;31m' 
set Purple='\033[0;35m'
set Color_Off='\033[0m'


if( $#argv != 6) then
echo "${Red}   Run_it_all.csh num                                          ${Color_Off}"
echo "${Red}   num       = job number                                      ${Color_Off}"
echo "${Red}   type of sim = single, pisa, embed, all sim (no embed), all  ${Color_Off}"
echo "${Red}   type of input  = single, helios, pythia8                    ${Color_Off}"
echo "${Red}   type of particle  = photon, other                           ${Color_Off}"
echo "${Red}   shift  = 1-10000 - Nev                                      ${Color_Off}"
echo "${Red}   Number of events  = 1-10000                                 ${Color_Off}"
exit -1
endif

set jobno = $1
set selected_paticle = $4
set shift = $5
set NEVT = $6

echo "${Green} jobno is set to     $jobno              ${Color_Off}"
echo "${Green} type of sim         $2                  ${Color_Off}"
echo "${Green} type of input       $3                  ${Color_Off}"
echo "${Green} selected_paticle    $selected_paticle   ${Color_Off}"
echo "${Green} shift               $shift              ${Color_Off}"
echo "${Green} NEVT                $NEVT               ${Color_Off}"

setenv DATADIR $PWD

source /opt/phenix/bin/phenix_setup.csh new
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH $DATADIR/embedreco/install/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH $DATADIR/embed/svx_cent_ana/install/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH .:/phenix/hhj/hachiya/15.08/embed/embed:$TSEARCHPATH
setenv TSEARCHPATH .:/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/install/:$TSEARCHPATH

echo "${Green} LD_LIBRARY_PATH: ${Color_Off}"
echo $LD_LIBRARY_PATH
echo "${Green} TSEARCHPATH: ${Color_Off}"
echo $TSEARCHPATH

set outputsingle_dir = $DATADIR/output_single

set INPUT = `expr $shift + $jobno`
#echo "${Green} input is " $INPUT
set DIR = `printf "%05d" $INPUT`

echo "${Green} INPUT               $INPUT              ${Color_Off}"
echo "${Green} DIR                 $DIR                ${Color_Off}"
echo "${Green} outputsingle_dir    $outputsingle_dir   ${Color_Off}"

if( -d $outputsingle_dir ) then
echo "${Green} outputsingle_dir exists ${Color_Off}"
else 
echo "${Green} mkdir $outputsingle_dir ${Color_Off}"
mkdir $outputsingle_dir
chmod g+rx $outputsingle_dir
endif

set RUNNUMBERS_RG0 = (409471 409471 409471)
#set RUNNUMBERS_RG0 = (409149 409150 409147 409151 409086 409087 409088 409089 409092 409115 409116 409118 409120 409121 409123 409124 409125 409152 409153 409154 409200 409222 409223 409224 409225 409226 409227 409228 409229 409230 409231 409232 409233 409235 409236 409237 409239 409241 409242 409243 409244 409295 409297 409298 409299 409300 409301 409302 409303 409304 409311 409312 409313 409305 409306 409307 409308 409310 409436 409437 409438 409428 409429 409430 409431 409432 409433 409435 409439 409440 409442 409443 409454 409455 409456 409457 409458 409459 409465 409467 409469 409471)
set RANDOM=$$
set rnd = `expr $RANDOM % ${#RUNNUMBERS_RG0} + 1` ### for bash use like this and numering is from 0 ${#RUNNUMBERS_RG0[@]}
set run_number = $RUNNUMBERS_RG0[$rnd]

set input_real_dir = $DATADIR/real/work/outputfull
set runlistname = $outputsingle_dir/runs_$jobno.txt
ls $input_real_dir/vertexes_MB* > $runlistname
set list = `cat $runlistname`
set RUNNUMBERS = ( )
foreach file ( $list )

  set runnum = `basename $file .PRDFF| awk -F'-' '{printf "%d", $2}'`
  set seqnum = `basename $file .PRDFF| awk -F'-' '{printf "%s", $3}'`
  set RUNNUMBERS = ( $RUNNUMBERS $runnum )
end

@ irun = $jobno % $#RUNNUMBERS + 1
set run_number = $RUNNUMBERS[$irun]
set vertex_txt_dir = $DATADIR/real/work/outputfull/vertexes_MB-0000$run_number-0001.txt

echo "${Green} $irun $RUNNUMBERS $#RUNNUMBERS                 ${Color_Off}"
echo "${Green} =======Generating a Random Run Number==========${Color_Off}"
echo "${Green}                 $run_number                    ${Color_Off}"

echo "${Green} vertex_txt_dir  $vertex_txt_dir                ${Color_Off}"
rm $runlistname

echo "${Purple}===============================================${Color_Off}"
echo "${Purple}============= START SIMULATION  ===============${Color_Off}"
echo "${Purple}===============================================${Color_Off}"

set inputvtx = $vertex_txt_dir
set oscarname = $DIR.oscar.particles.dat

if( ( $2 == 0 || $2 == 3 || $2 == 4 ) && $3 == 0 ) then
echo "${Purple}===============================================${Color_Off}"
echo "${Purple}============= SINGLE SIM STARTS ===============${Color_Off}"
echo "${Purple}===============================================${Color_Off}"
set scriptdir   = $DATADIR/sim/gen/single
set macroname   = make_single.C
set outsingle   = $DATADIR/output_single/single
set tmpdir      = "/phenix/plhf/${USER}/tmp/job_single_$INPUT"
set ptmin = 0.4
set ptmax = 10.0
set n     = -1. #n: <0 hagdorn (mb HeAu), =0 flat, >0 power law
set id    = $selected_paticle #0,1,2,3,4,5,6-pi0,pi+,pi-,e-,e+,p,antip################helios jpsi and phi####pythia ccbar bbar

echo "${Green}jobno          $jobno            ${Color_Off}"
echo "${Green}run_number     $run_number       ${Color_Off}"
echo "${Green}ptmin          $ptmin            ${Color_Off}"
echo "${Green}ptmax          $ptmax            ${Color_Off}"
echo "${Green}n              $n                ${Color_Off}"
echo "${Green}id             $id               ${Color_Off}"
echo "${Green}inputvtx       $inputvtx         ${Color_Off}"
echo "${Green}oscarname      $oscarname        ${Color_Off}"
echo "${Green}scriptdir      $scriptdir        ${Color_Off}"
echo "${Green}macroname      $macroname        ${Color_Off}"
echo "${Green}outsingle      $outsingle        ${Color_Off}"
echo "${Green}tmpdir         $tmpdir           ${Color_Off}"

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "${Green}cd $tmpdir${Color_Off}"
cd       $tmpdir

cp $scriptdir/$macroname .
echo "${Green}start running root -l -b -q" 'make_single.C("'$vertex_txt_dir'","'$oscarname'",'$NEVT','$ptmin','$ptmax','$n','$id')'"${Color_Off}"
root -l -b -q 'make_single.C("'$vertex_txt_dir'","'$oscarname'",'10000','$ptmin','$ptmax','$n','$id')'
echo "${Green}finished running${Color_Off}"
echo "${Green}copying${Color_Off}"
cp $oscarname  $outsingle/
echo "${Green}finished copying${Color_Off}"
#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "${Green}removed $tmpdir${Color_Off}"

endif

if( ( $2 == 0 || $2 == 3 || $2 == 4 ) && $3 == 1 ) then
echo "${Purple}===============================================${Color_Off}"
echo "${Purple}============= HELIOS TO OSCAR =================${Color_Off}"
echo "${Purple}===============================================${Color_Off}"
set inputhelios = $DATADIR/output_single/helios/helios_phi_ee_02_5_10M.root
if( $selected_paticle == 2) then
set inputhelios = $DATADIR/output_single/helios/helios_jpsi_ee_0_10_50M.root
 #set inputhelios = $DATADIR/output_single/helios/helios_jpsi_ee_02_5_11M.root
endif
if( $selected_paticle == 0) then
 set inputhelios = $DATADIR/output_single/helios/helios_pi0_gg_04_6_100M.root
endif
if( $selected_paticle == 3) then
 set inputhelios = $DATADIR/output_single/helios/helios_pi0_gee_05_2_10M.root
endif
set scriptdir   = $DATADIR/sim/gen/HELIOS/work
set scriptname  = Convert_HELIOS.csh
set macroname   = WriteROOT2Oscar.C
set outsingle   = $DATADIR/output_single/helios
set tmpdir      = "/phenix/plhf/${USER}/tmp/job_helios_$INPUT"

echo "${Green}jobno          $jobno            ${Color_Off}"
echo "${Green}run_number     $run_number       ${Color_Off}"
echo "${Green}inputhelios    $inputhelios      ${Color_Off}"
echo "${Green}inputvtx       $inputvtx         ${Color_Off}"
echo "${Green}oscarname      $oscarname        ${Color_Off}"
echo "${Green}scriptdir      $scriptdir        ${Color_Off}"
echo "${Green}scriptname     $scriptname       ${Color_Off}"
echo "${Green}macroname      $macroname        ${Color_Off}"
echo "${Green}outsingle      $outsingle        ${Color_Off}"
echo "${Green}tmpdir         $tmpdir           ${Color_Off}"

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "${Green}cd $tmpdir${Color_Off}"
cd       $tmpdir

cp $scriptdir/$scriptname .
cp $scriptdir/$macroname .
echo "${Green}running ./$scriptname $inputvtx $jobno 10000 $inputhelios $oscarname${Color_Off}"
./$scriptname $inputvtx $jobno 10000 $inputhelios $oscarname
echo "${Green}finished running${Color_Off}"
echo "${Green}copying${Color_Off}"
cp $oscarname  $outsingle/
echo "${Green}finished copying${Color_Off}"
#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "${Green}removed $tmpdir${Color_Off}"

endif

if( ( $2 == 0 || $2 == 3 || $2 == 4 ) && $3 == 2 ) then
echo "${Purple}===============================================${Color_Off}"
echo "${Purple}============= PYTHIA TO OSCAR =================${Color_Off}"
echo "${Purple}===============================================${Color_Off}"
set inputpythia = $DATADIR/output_single/pythia8/ccbartree$DIR.root
if( $selected_paticle == 2) then
 set inputpythia = $DATADIR/output_single/pythia8/bbbartree$DIR.root
endif
if( $selected_paticle == 3) then
 set inputpythia = $DATADIR/output_single/pythia8/jetpairstree$DIR.root
endif
set scriptdir   = $DATADIR/sim/gen/pythia8
set scriptname  = Convert_pythia8.csh
set macroname   = WriteROOT2OscarPythia.C
set outsingle   = $DATADIR/output_single/pythia8
set tmpdir      = "/phenix/plhf/${USER}/tmp/job_pythia_$INPUT"

echo "${Green}jobno          $jobno            ${Color_Off}"
echo "${Green}run_number     $run_number       ${Color_Off}"
echo "${Green}inputpythia    $inputpythia      ${Color_Off}"
echo "${Green}inputvtx       $inputvtx         ${Color_Off}"
echo "${Green}oscarname      $oscarname        ${Color_Off}"
echo "${Green}scriptdir      $scriptdir        ${Color_Off}"
echo "${Green}scriptname     $scriptname       ${Color_Off}"
echo "${Green}macroname      $macroname        ${Color_Off}"
echo "${Green}outsingle      $outsingle        ${Color_Off}"
echo "${Green}tmpdir         $tmpdir           ${Color_Off}"

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "${Green}cd $tmpdir${Color_Off}"
cd       $tmpdir

cp $scriptdir/$scriptname .
cp $scriptdir/$macroname .
echo "${Green}running ./$scriptname $inputvtx $inputpythia $oscarname ${Color_Off}"
./$scriptname $inputvtx $inputpythia $oscarname
echo "${Green}finished running${Color_Off}"
echo "${Green}copying${Color_Off}"
cp $oscarname  $outsingle/
echo "${Green}finished copying${Color_Off}"
#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "${Green}removed $tmpdir${Color_Off}"

endif


if( $2 == 1 || $2 == 3 || $2 == 4 ) then

echo "${Purple}===============================================${Color_Off}"
echo "${Purple}============= START PISA NOW ==================${Color_Off}"
echo "${Purple}===============================================${Color_Off}"

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
set tmpdir      = "/phenix/plhf/${USER}/tmp/job_pisa_$INPUT"

echo "${Green}jobno          $jobno            ${Color_Off}"
echo "${Green}run_number     $run_number       ${Color_Off}"
echo "${Green}Nev            $NEVT             ${Color_Off}"
echo "${Green}oscarname      $oscarname        ${Color_Off}"
echo "${Green}pisadir        $pisadir          ${Color_Off}"
echo "${Green}outpisadir     $outpisadir       ${Color_Off}"
echo "${Green}tmpdir         $tmpdir           ${Color_Off}"

#checkng outdir
if( ! -d $outpisadir ) then
echo "${Green}creating $outpisadir${Color_Off}" 
mkdir -p $outpisadir
endif
#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "${Green}cd $tmpdir${Color_Off}"
cd       $tmpdir


cp $oscarname oscar.particles.dat
cp $pisadir/* .

sed -i 's/ptrig [0-9]*/ptrig '$NEVT'/' glogon.kumac

echo "${Green}running pisa<pisa.input${Color_Off}"
pisa<pisa.input
echo "${Green}finished running${Color_Off}"

echo "${Purple}===============================================${Color_Off}"
echo "${Purple}================ PISA TO DST ==================${Color_Off}"
echo "${Purple}===============================================${Color_Off}"

echo "${Green}running root -b -q -l pisaToDST_VTX.C${Color_Off}"
root -b -q -l 'pisaToDST_VTX.C('$NEVT', "PISAEvent.root", "dst_out.root", "svxeval.root", '$run_number', 0)'
echo "${Green}finished running${Color_Off}"

if ( -f dst_out.root ) then
echo "${Green}moving${Color_Off}"
mv dst_out.root $outpisadir/$outdstname
echo "${Green}finished moving${Color_Off}"
else 
echo "${Green}something is fkng wrong${Color_Off}"
endif

#remove tmp dir
cd $DATADIR
rm -fr $tmpdir
echo "${Green}removed $tmpdir${Color_Off}"

endif

##########################################
################EMBEDDING#################
#################STARTS###################
##########################################
if( $2 == 2 || $2 == 4 ) then
echo "${Purple}===============================================${Color_Off}"
echo "${Purple}================= EMBEDDING ===================${Color_Off}"
echo "${Purple}===============================================${Color_Off}"

echo "${Green}Compiling libs${Color_Off}"
setenv sourcedir /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/embed/svx_cent_ana/build
if( $1 == 0) then
cd $sourcedir
make -j4
make install
endif

echo "${Green}End of Compiling${Color_Off}"
echo "${Green}Setting internal variables...${Color_Off}"

set runnum    = $run_number
set jobno     = $1
set evtnum    = $NEVT

set INPUT = `expr $shift + $jobno`
#echo "${Green}$INPUT ${Color_Off}"
set outntana  = tree"$DIR".root

if ( $3 == 1 ) then 
set inputoscar = $DATADIR/output_single/helios/$DIR.oscar.particles.dat
else if ( $3 == 2 ) then 
set inputoscar = $DATADIR/output_single/pythia8/$DIR.oscar.particles.dat
else 
set inputoscar = $DATADIR/output_single/single/$DIR.oscar.particles.dat
endif
set inputsim = $DATADIR/output_single/simdst/dst_out_single_$DIR.root
#set inputvtx = $DATADIR/real/work/output/vertexes.txt
set inputreal = $DATADIR/real/work/outputfull/CNTmerge_MB-0000$run_number-0001.root
set outdst    = kek0.root
set embedpartID = $selected_paticle

set scriptdir = $DATADIR/embed/work

set outmytreedir   = $DATADIR/output_single/embed

#checkng outdir
if( ! -d $outmytreedir ) then
echo "${Green}creating $outmytreedir${Color_Off}" 
mkdir -p $outmytreedir
endif

set tmpdir      = "/phenix/plhf/${USER}/tmp/job_embed_$INPUT"

set inreal = `basename $inputreal`
set insim  = `basename $inputsim`

set runnum = `echo $inreal | awk -F'-' '{printf "%d", $2}'`

echo "${Green}runnum       $runnum         ${Color_Off}"
echo "${Green}jobno        $jobno          ${Color_Off}"
echo "${Green}evtnum       $evtnum         ${Color_Off}"
echo "${Green}inputreal    $inputreal      ${Color_Off}"
echo "${Green}inputsim     $inputsim       ${Color_Off}"
echo "${Green}outdst       $outdst         ${Color_Off}"
echo "${Green}embedpartID  $embedpartID    ${Color_Off}"
echo "${Green}outntana     $outntana       ${Color_Off}"
echo "${Green}inputvtx     $inputvtx       ${Color_Off}"
echo "${Green}inputoscar   $inputoscar     ${Color_Off}"
echo "${Green}scriptdir    $scriptdir      ${Color_Off}"
echo "${Green}outmytreedir $outmytreedir   ${Color_Off}"
echo "${Green}tmpdir       $tmpdir         ${Color_Off}"

# exit


#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "${Green}cd $tmpdir${Color_Off}"
cd       $tmpdir


##pisaToDSTLinker with new library
/opt/phenix/core/bin/LuxorLinker.pl -1 $runnum

cp ${scriptdir}/Fun4All_embedeval.C .
cp ${scriptdir}/Fun4All_embedeval_svx.C .
cp ${scriptdir}/embed_IOManager.C    .
cp ${scriptdir}/svxPISA.par         .
cp ${scriptdir}/svx_threshold.dat   .

#copy input file
#cp $inputreal .
#cp $inputsim .
#echo "pwd"
#pwd
#ls -ltr
#echo "yolo"
echo "${Green}running in root: .x Fun4All_embedeval_svx.C($evtnum, "'"'$insim'"'", "'"'$inreal'"'", "'"'$outdst'"'", "$embedpartID", "'"'$outntana'"'", "'"'$inputvtx'"'", "'"'$inputoscar'"'", $runnum);${Color_Off}"
echo ".x Fun4All_embedeval_svx.C($evtnum, "'"'$inputsim'"'", "'"'$inputreal'"'", "'"'$outdst'"'", "$embedpartID", "'"'$outntana'"'", "'"'$inputvtx'"'", "'"'$inputoscar'"'", $runnum);" >  cmd.input
echo ".q" >> cmd.input

##run root
root -b < cmd.input #>& $HOME/root.log
echo "${Green}end of running${Color_Off}"
#ls -ltr
#move
#mv -f $outdst    $outdstdir
#mv -f kek2.root $outmytreedir/
#mv -f $outntana  $outmytreedir/

echo "${Green}start moving outfiles${Color_Off}"
mv -f my-$outntana $outmytreedir/
echo "${Green}end of moving outfiles${Color_Off}"

#remove tmp dir
echo "${Green}start removing tmp dir${Color_Off}"
cd $HOME
rm -fr $tmpdir
echo "${Green}removed: $tmpdir${Color_Off}"

endif

echo "${Purple}===============================================${Color_Off}"
echo "${Purple}================== THE END ====================${Color_Off}"
echo "${Purple}===============================================${Color_Off}"
##########################################
###################END####################
################EMBEDDING#################
##########################################