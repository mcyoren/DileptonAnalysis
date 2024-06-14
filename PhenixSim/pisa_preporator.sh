RUNNO=406541
run=14

mkdir make_sim
mkdir make_sim/pisa

cd make_sim/pisa

LuxorLinker.pl -1 $RUNNO
/afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run14/pisaLinker.csh
/afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run14/pisaToDSTLinker.csh
./pisaToDSTLinker.csh

cp /afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run14/pisaToDST_200GeV.C pisaToDST.C

cp ../../pisa_svx/pisaToDST_VTX.C .
cp ../../pisa_svx/glogon.kumac .
cp ../../pisa_svx/pisa.input .
cp ../../pisa_svx/pisa.kumac .  ## use your un number!!!!!!
cp ../../pisa_svx/pisaToDST_IOManager.C .
cp ../../pisa_svx/svx* .

#cp ../glogon.kumac .
#cp ../pisa.input .
#cp ../pisaToDST.C .
