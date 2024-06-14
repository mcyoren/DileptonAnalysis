
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs/:/phenix/u/vdoomra/install
export INSTALL=/phenix/u/vdoomra/install
export LD_LIBRARY_PATH=$INSTALL/lib:$LD_LIBRARY_PATH
export TSEARCHPATH=/phenix/u/vdoomra/install

cd output/00000/pisa/

#cp /phenix/plhf/vdoomra/simulation/run_TTreeMaker.C .

root -l -b -q 'run_TTreeMaker.C("dst_out.root", "output.root", 0, 0)'

