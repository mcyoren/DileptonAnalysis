#!/bin/csh

source /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/etc/eic_cshrc.csh -n
setenv PYTHIA8DATA /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/share/Pythia8/xmldoc

make WriteTreeoutputccbar
make main113_pA_nPDF

if( $#argv != 0) then
 if( $1 == 0) then
  ./WriteTreeoutputccbar
 endif
 if( $1 == 1) then
  python3 main113_pA_nPDF.py
 endif
endif
