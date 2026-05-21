#!/bin/csh

setenv EIC_LEVEL EIC2022a
source /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/etc/eic_cshrc.csh

set shift = $4
set INPUT = `expr $shift + $2`
echo $INPUT

set tmpdir = "/phenix/plhf/mitran/tmp/${USER}_job_pythia_$INPUT"
set sourcedir = /phenix/plhf/mitran/Simul/Dileptons/sim/gen/pythia8
set outputdir = /phenix/plhf/mitran/Simul/Dileptons/output_single/pythia8
set DIR = `printf "%05d" $INPUT`
set executable = WriteTreeoutputccbar
set outname = ccbar


if( $1 == 1) then
 set executable = WriteTreeoutputbbbar
 set outname = bbbar
endif

if( $1 == 2) then
 set executable = WriteTreeOutputJetPairs
 set outname = jetpairs
endif

if( $1 == 4) then
 set executable = WriteTreeoutputccbarSoft
 set outname = ccbarSoft1
endif

if( $1 == 5) then
 set executable = WriteLambda
 set outname = Lambda
endif

if( $1 == 6) then
 set executable = WriteTreeoutputccbarSoft
 set outname = ccbarLambda
endif

if( $1 == 7) then
 set executable = WriteTreeoutputccbarSoft
 set outname = ccbarLambdaDetroit
endif

 
echo "all params are set to"

echo "tmpdir               $tmpdir    "
echo "job                  $2         "
echo "sourcedir            $sourcedir "
echo "outputdir            $outputdir "
echo "DIR                  $DIR       "
echo "seed                 $INPUT     "
echo "Nev                  $3         "
echo "shift                $4         "
echo "running executable:  $executable"
echo "outname:             $outname   "

if( $2 == 0) then
 cd $sourcedir
 make $executable
endif

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir
cp $sourcedir/$executable* .
cp $sourcedir/main113_pA_nPDF* .
cp $sourcedir/Makefile* .

if( $#argv != 0) then
 if( $1 == 0) then
  ./$executable $INPUT $3
 endif
 if( $1 == 1) then
  ./$executable $INPUT $3
 endif
 if( $1 == 2) then
  ./$executable $INPUT $3
 endif
 if( $1 == 3) then
  make main113_pA_nPDF
  python3 main113_pA_nPDF.py
 endif
 if( $1 == 4) then
  ./$executable $INPUT $3
 endif
 if( $1 == 5) then
  ./$executable $INPUT $3
 endif
 if( $1 == 6) then
  ./$executable $INPUT $3
 endif
 if( $1 == 7) then
  ./$executable $INPUT $3
 endif

endif

echo "cp *.root $outputdir/""$outname""newtree$DIR.root"
mv *.root $outputdir/"$outname"newtree$DIR.root
#cp *.dat  $outputdir/$DIR.oscar.parcticles.dat
#remove tmp dir
cd $HOME
rm -fr $tmpdir
echo "removed $tmpdir"
