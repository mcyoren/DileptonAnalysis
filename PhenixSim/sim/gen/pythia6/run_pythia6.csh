#!/bin/csh

source /opt/phenix/core/bin/phenix_setup.csh -n

setenv PYTHIA6_DIR /cvmfs/phenix.sdcc.bnl.gov/x8664_sl7/opt/phenix/core/stow/pythia-6.4.28
setenv LD_LIBRARY_PATH ${PYTHIA6_DIR}/lib:${LD_LIBRARY_PATH}

gfortran -m32 -O3 -o run_pythia6_mb_ccbar_dielectrons_forced \
  run_pythia6_mb_ccbar_dielectrons_forced.f \
  -L${PYTHIA6_DIR}/lib -lPythia6

set JOBID = $1
set TARGET_PAIRS = $2

./run_pythia6_mb_ccbar_dielectrons_forced ${TARGET_PAIRS} ${JOBID}