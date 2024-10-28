#!/bin/csh

echo "hello"

source /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/etc/eic_cshrc.csh -n
setenv PYTHIA8DATA /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/share/Pythia8/xmldoc

# kek.root 500 phi phi-\>ee 0.2 10.
echo "name:kek.root Nev:500 part:phi decay:phi-\>ee pt_min:0.2 pt_max:10."


root -l -b << EOF
    .L ../source/HELIOSLibrary/WriteEvent.C+
    .L WriteROOTOutput.C+
EOF

if( $#argv != 0) then
 if( $1 == 0) then
 echo "running all standard params"
  root -l -b << EOF
    gSystem->Load("../source/HELIOSLibrary/WriteEvent_C.so")
    .x WriteROOTOutput.C+
EOF
 endif
 if( $1 != 0) then
 
 set outfile   = $1
 set nev       = $2
 set partname  = $3
 set decayname = $4
 set ptmin     = $5
 set ptmax     = $6
 
 echo "all params are set to"
 
 echo "outfile      $outfile    "
 echo "nev          $nev        "
 echo "partname     $partname   "
 echo "decayname    $decayname  "
 echo "ptmin        $ptmin      "
 echo "ptmax        $ptmax      "
 echo "running WriteROOTOutput("$outfile",$nev,"$partname","$decayname",$ptmin,$ptmax)"

  root -l -b << EOF
    gSystem->Load("../source/HELIOSLibrary/WriteEvent_C.so")
    gSystem->Load("WriteROOTOutput_C.so")
    WriteROOTOutput("$outfile",$nev,"$partname","$decayname",$ptmin,$ptmax)
EOF
 endif
endif

echo "the end!"
