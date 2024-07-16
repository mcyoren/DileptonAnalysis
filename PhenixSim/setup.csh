#!/bin/csh

if( $#argv != 1) then
echo "reco_rawhit.csh num"
echo "   num       = job number"
echo "   evtnum    = Nevent"
echo "   inputreal = inputfile"
echo "   inputsim  = inputfile"
echo "   outdst    = outdst"
echo "   outntuple = outdst"
echo "   outntana  = outdst"
exit -1
endif

if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir

cd sim
./pida_preporator.csh

./run_sim_single.csh

condor_submit runs_sim_single.job

cd $PWD/real
