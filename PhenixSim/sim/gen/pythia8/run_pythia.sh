#!/bin/sh

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR

cp /phenix/plhf/vdoomra/standalone_pythia8/WriteTreeoutputbbbar .

./WriteTreeoutputbbbar

cp *.root /phenix/plhf/vdoomra/standalone_pythia8/output_bottom_single_electron/bottom_$1.root

popd
rm -r $DIR
