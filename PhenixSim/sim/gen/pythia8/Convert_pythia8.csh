#!/bin/csh

source /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/etc/eic_cshrc.csh -n
setenv PYTHIA8DATA /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/share/Pythia8/xmldoc

root -l -b << EOF
 .x WriteROOT2OscarPythia.C("$1","$2","$3")
EOF