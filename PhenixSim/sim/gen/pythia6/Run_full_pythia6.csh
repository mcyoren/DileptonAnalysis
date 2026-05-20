#!/bin/csh
# source the alma9 setup script to get our alma9 environment
source /opt/phenix/core/bin/phenix_alma9_setup.csh -n
# run our SL7 script in an SL7 singularity container

echo "STARTING Run_full_pythia6.csh"
set sart_time = `date`
echo "Running on alma9, but using SL7 singularity container for pythia6 and TTreeMaker"

set Nev = $1
set shift = $3
set INPUT = `expr $shift + $2`
echo $INPUT

set tmpdir = "/home/tmp/${USER}_job_pythia6_$INPUT"
set sourcedir = /phenix/plhf/mitran/Simul/Dileptons/sim/gen/pythia6
set outputdir = /phenix/plhf/mitran/Simul/Dileptons/output_single/pythia6
set DIR = `printf "%05d" $INPUT`
set pythia6_macro = run_pythia6_mb_ccbar_dielectrons_forced.f
set TTreeMaker_macro = ConvertPythia6PairsToTree.cc
set pythia6_script = Run_pythia6_alma9.csh
set ttreemaker_script = Run_TTreeMaker_alma9.csh
set pythia6_script_sl7 = run_pythia6.csh
set ttreemaker_script_sl7 = run_ttree_maker.csh
set outname = ccbar
 
echo "all params are set to"

echo "tmpdir           $tmpdir           "
echo "job              $2                "
echo "sourcedir        $sourcedir        "
echo "outputdir        $outputdir        "
echo "DIR              $DIR              "
echo "seed             $INPUT            "
echo "Nev              $3                "
echo "shift            $4                "
echo "running:         $pythia6_script   "
echo "running:         $ttreemaker_script"
echo "outname:         $outname          "

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir
cp $sourcedir/$pythia6_macro .
cp $sourcedir/$TTreeMaker_macro .
cp $sourcedir/$pythia6_script .
cp $sourcedir/$ttreemaker_script .
cp $sourcedir/$pythia6_script_sl7 .
cp $sourcedir/$ttreemaker_script_sl7 .

echo "running pythia6 script"
./$pythia6_script $INPUT $Nev
echo "Done running pythia6 script"
echo "running TTreeMaker script"
./$ttreemaker_script "$outname"newtree$DIR.root
echo "Done running TTreeMaker script"

echo "mv *.root $outputdir/""$outname""newtree$DIR.root"
mv *.root $outputdir/"$outname"newtree$DIR.root
#cp *.dat  $outputdir/$DIR.oscar.parcticles.dat
#remove tmp dir
cd $HOME
rm -fr $tmpdir
echo "removed $tmpdir"

set end_time = `date`
set start_sec = `date -d "$sart_time" +%s`
set end_sec = `date -d "$end_time" +%s`
set duration_sec = `expr $end_sec - $start_sec`
set duration_min = `expr $duration_sec / 60`
set duration_hr = `expr $duration_min / 60`
echo "Job duration: $duration_sec seconds, or $duration_min minutes, or $duration_hr hours"
echo "FINISHED Run_full_pythia6.csh"

