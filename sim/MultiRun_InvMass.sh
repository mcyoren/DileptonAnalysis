#!/bin/bash

Purple='\033[0;35m'
NC='\033[0m' # No Color

if [ "$#" -ne 5 ]; then
    echo -e "${Purple}You need to run the script with the following 4 parameters: input file, number of threads, output file name, number of events${NC}, particle type : pi0, phi, jpsi, ccbar, bbbar, omega, qgp, psi2s"
    exit 1
fi

echo -e "${Purple}input file: $1\nnumber of threads: $2\noutput file name: $3\nnumber of events: $4${NC}"

root -l -b << EOF
    .L MyEvent_OLD.C+
    .L ../AnaTrain/Run14AuAuLeptonComby/MyML.C+
    .L ../AnaTrain/Run14AuAuLeptonComby/MyEvent.C+
    .L InvMass.C+
EOF

echo -e "${Purple}compiling is done, starting running${NC}"

for i in `seq 1 $2`
do
    ./run_InvMass.sh $1 $i $2 $4 $5 &
done
wait

echo -e "${Purple}running is done, starting merging${NC}"

hadd -f -k -O $3 my-kek*
rm my-kek*
echo -e "${Purple}All done!${NC}"