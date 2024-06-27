# PYTHIA configuration file.
# Generated on Tue Jun 28 11:41:02 EDT 2022 with the user supplied options:
# --prefix=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a
# --with-fastjet3
# --with-lhapdf6
# --with-lhapdf6-lib=/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/LHAPDF-6.2.3//lib
# --with-lhapdf6-include=/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/LHAPDF-6.2.3//include
# --with-hepmc3
# --with-gzip

# Install directory prefixes.
PREFIX_BIN=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/bin
PREFIX_INCLUDE=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/include
PREFIX_LIB=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/lib
PREFIX_SHARE=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/gcc/8.3.0.1-0a5ad/x86_64-centos7/bin/g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -DGZIP
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=
FASTJET3_INCLUDE=-I/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/fastjet-3.3.2/include
FASTJET3_LIB=-L/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/CGAL-4.14.0/lib -Wl,-rpath,/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/CGAL-4.14.0/lib -L/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/boost-1.70.0/lib -Wl,-rpath,/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/boost-1.70.0/lib -L/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/fastjet-3.3.2/lib -Wl,-rpath,/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/fastjet-3.3.2/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=
HEPMC3_INCLUDE=-I/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/release_new/new.4/include
HEPMC3_LIB=-L/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/release_new/new.4/lib64 -Wl,-rpath,/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/release_new/new.4/lib64 -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=
LHAPDF6_INCLUDE=-I/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/LHAPDF-6.2.3/include
LHAPDF6_LIB=-L/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/LHAPDF-6.2.3/lib -Wl,-rpath,/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/LHAPDF-6.2.3/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=true
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=-lz

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
