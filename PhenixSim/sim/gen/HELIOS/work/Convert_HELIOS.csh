#!/bin/csh

echo "hello"

source /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/etc/eic_cshrc.csh -n
setenv PYTHIA8DATA /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/share/Pythia8/xmldoc


root -l -b << EOF
    .L /phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/WriteEvent.C+
    .L WriteROOT2Oscar.C+
EOF

echo "lib is builded"

if( $#argv != 0) then
 if( $1 == 0) then
 echo "running all standard params"
  root -l -b << EOF
    gSystem->Load("/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/WriteEvent_C.so")
    gSystem->Load("WriteROOT2Oscar_C.so")
    WriteROOT2Oscar()
EOF
 endif
 if( $1 != 0) then
 
 set vtxfile   = $1
 set iter      = $2
 set nev       = $3
 set infile    = $4
 set outfile   = $5
 
 echo "all params are set to"
 
 echo "vtxfile    $vtxfile  "
 echo "iter       $iter     "
 echo "nev        $nev      "
 echo "infile     $infile   "
 echo "outfile    $outfile  "
 echo "running WriteROOT2Oscar("$vtxfile","$iter",$nev,"$infile","$outfile")"

  root -l -b << EOF
    gSystem->Load("/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/WriteEvent_C.so")
    gSystem->Load("WriteROOT2Oscar_C.so")
    WriteROOT2Oscar("$vtxfile",$iter,$nev,"$infile","$outfile")
EOF
 endif
endif

echo "the end!"
