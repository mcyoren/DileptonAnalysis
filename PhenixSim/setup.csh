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
